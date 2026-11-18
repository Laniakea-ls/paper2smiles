"""Persistent, user-isolated PDF jobs behind the Portal signed-identity gateway."""
import asyncio,base64,hashlib,hmac,json,os,re,signal,sqlite3,sys,time,uuid
from contextlib import asynccontextmanager
from pathlib import Path
from urllib.parse import unquote
from starlette.applications import Starlette
from starlette.requests import Request
from starlette.responses import JSONResponse,FileResponse
from starlette.routing import Route
from starlette.exceptions import HTTPException
from paper2smiles.settings import ROOT,settings

STAGES={'plan':'识别文献','activity_ocr':'提取活性表','si_discovery':'查找补充材料','names':'提取完整名称','opsin':'名称转结构','reference_names':'核对名称','markush':'通式与片段连接','full_structure':'结构图识别','report':'生成报告','completed':'完成'}

def root():return Path(settings()['jobs_root']).resolve()
def connect():
    root().mkdir(parents=True,exist_ok=True)
    con=sqlite3.connect(root()/'jobs.sqlite3',timeout=20);con.row_factory=sqlite3.Row
    con.execute('PRAGMA journal_mode=WAL')
    con.execute('CREATE TABLE IF NOT EXISTS jobs (id TEXT PRIMARY KEY,owner TEXT NOT NULL,filename TEXT NOT NULL,pdf_sha TEXT NOT NULL,status TEXT NOT NULL,created REAL NOT NULL,updated REAL NOT NULL,error TEXT NOT NULL DEFAULT \'\')')
    return con

def owner(request):
    token=request.headers.get('authorization','').removeprefix('Bearer ')
    secret=settings()['identity_secret']
    try:
        header,payload,signature=token.split('.')
        decode=lambda s:base64.urlsafe_b64decode(s+'='*(-len(s)%4))
        if len(secret)<32:raise ValueError()
        if json.loads(decode(header)).get('alg')!='HS256':raise ValueError()
        expected=hmac.new(secret.encode(),(header+'.'+payload).encode(),hashlib.sha256).digest()
        if not hmac.compare_digest(expected,decode(signature)):raise ValueError()
        claims=json.loads(decode(payload));now=time.time()
        if claims.get('iss')!='aidd-portal' or claims.get('aud')!=settings()['audience'] or claims.get('exp',0)<=now or claims.get('iat',now+100)>now+10:raise ValueError()
        if claims.get('path')!=request.url.path or claims.get('method')!=request.method or not isinstance(claims.get('sub'),str) or not claims['sub']:raise ValueError()
        return claims['sub']
    except (ValueError,TypeError,KeyError,json.JSONDecodeError):raise HTTPException(401,'Authentication required')

def job_for(request):
    user=owner(request);jid=request.path_params['job_id']
    if not re.fullmatch('[a-f0-9]{32}',jid):raise HTTPException(404,'任务不存在')
    with connect() as con:row=con.execute('SELECT * FROM jobs WHERE id=? AND owner=?',(jid,user)).fetchone()
    if not row:raise HTTPException(404,'任务不存在')
    return dict(row)

def public_job(row):
    row=dict(row);jid=row['id'];folder=root()/jid
    value={k:row[k] for k in ('id','filename','status','created','updated','error')}
    status=folder/'status.json'
    if status.exists():
        try:
            status=json.loads(status.read_text());value['stage']=status.get('stage','');value['message']=STAGES.get(status.get('stage'),str(status.get('message',''))[:160])
        except (OSError,ValueError):pass
    summary=folder/'summary.json'
    if summary.exists():
        try:
            s=json.loads(summary.read_text());value['summary']={k:s[k] for k in ['title','doi','compounds','with_structure','measurements'] if k in s}
        except (OSError,ValueError):pass
    if row['status']=='completed':value['report_path']=f'/api/paper-sar/jobs/{jid}/files/report-data.json'
    return value

async def health(request):
    owner(request);cfg=settings()
    return JSONResponse({'ok':True,'name':'文献活性提取','configured':bool(cfg['api_key']),'version':'1.0.0','max_upload_bytes':cfg['max_upload_bytes']})

async def jobs(request):
    user=owner(request)
    if request.method=='GET':
        with connect() as con:rows=con.execute('SELECT * FROM jobs WHERE owner=? ORDER BY created DESC LIMIT 100',(user,)).fetchall()
        return JSONResponse({'jobs':[public_job(r) for r in rows]})
    cfg=settings()
    if not cfg['api_key']:raise HTTPException(503,'模型服务尚未配置')
    filename=Path(unquote(request.query_params.get('filename','paper.pdf'))).name[:240]
    if not filename.lower().endswith('.pdf'):raise HTTPException(400,'请上传 PDF 文件')
    jid=uuid.uuid4().hex;folder=root()/jid;folder.mkdir(parents=True);pdf=folder/'upload.pdf';size=0;digest=hashlib.sha256()
    try:
        with pdf.open('wb') as stream:
            async for chunk in request.stream():
                size+=len(chunk)
                if size>cfg['max_upload_bytes']:raise HTTPException(413,'PDF 不能超过 50 MB')
                digest.update(chunk);stream.write(chunk)
        with pdf.open('rb') as stream: signature=stream.read(5)
        if size<5 or signature!=b'%PDF-':raise HTTPException(400,'文件不是有效 PDF')
        # Load XML consumers before PyMuPDF to avoid binary library conflicts.
        import openpyxl
        import fitz
        try:
            with fitz.open(pdf) as document:
                if document.needs_pass:raise HTTPException(400,'暂不支持加密 PDF')
                if not 0<len(document)<=cfg['max_pages']:raise HTTPException(400,'PDF 页数超出限制')
        except (RuntimeError,ValueError):raise HTTPException(400,'无法读取 PDF')
        with connect() as con:
            con.execute('BEGIN IMMEDIATE')
            existing=con.execute("SELECT * FROM jobs WHERE owner=? AND pdf_sha=? AND status IN ('queued','running','completed') ORDER BY created DESC LIMIT 1",(user,digest.hexdigest())).fetchone()
            if existing:
                pdf.unlink();folder.rmdir();return JSONResponse(public_job(existing))
            now=time.time();con.execute('INSERT INTO jobs VALUES (?,?,?,?,?,?,?,?)',(jid,user,filename,digest.hexdigest(),'queued',now,now,''))
        return JSONResponse(public_job({'id':jid,'filename':filename,'status':'queued','created':now,'updated':now,'error':''}),202)
    except BaseException:
        if pdf.exists():pdf.unlink()
        if folder.exists() and not any(folder.iterdir()):folder.rmdir()
        raise

async def detail(request):return JSONResponse(public_job(job_for(request)))

async def cancel(request):
    row=job_for(request)
    with connect() as con:con.execute("UPDATE jobs SET status='cancelled',updated=? WHERE id=? AND status IN ('queued','running')",(time.time(),row['id']))
    return JSONResponse({'ok':True})

async def retry(request):
    row=job_for(request)
    if row['status'] not in ('failed','cancelled'):raise HTTPException(409,'当前任务不能重试')
    with connect() as con:
        con.execute('BEGIN IMMEDIATE')
        if (root()/row['id']/'.worker-running').exists():raise HTTPException(409,'任务正在停止，请稍后重试')
        con.execute("UPDATE jobs SET status='queued',error='',updated=? WHERE id=? AND status IN ('failed','cancelled')",(time.time(),row['id']))
    return JSONResponse({'ok':True})

async def artifact(request):
    row=job_for(request)
    if row['status']!='completed':raise HTTPException(409,'报告尚未完成')
    base=(root()/row['id']/'report').resolve();relative=request.path_params['path'];path=(base/relative).resolve()
    allowed=relative=='report-data.json' or relative.startswith(('assets/','downloads/'))
    if not allowed or not path.is_relative_to(base) or not path.is_file() or path.suffix.lower() not in {'.json','.pdf','.xlsx','.csv','.tsv','.png','.svg','.jpg'}:raise HTTPException(404,'文件不存在')
    headers={'Cache-Control':'private, no-store','X-Content-Type-Options':'nosniff'}
    filename=path.name if path.suffix.lower() in {'.xlsx','.csv','.tsv'} else None
    return FileResponse(path,filename=filename,headers=headers)

async def worker_loop():
    for marker in root().glob('*/.worker-running'):marker.unlink()
    with connect() as con:con.execute("UPDATE jobs SET status='queued' WHERE status='running'")
    process=None
    try:
        while True:
            with connect() as con:
                con.execute('BEGIN IMMEDIATE');row=con.execute("SELECT * FROM jobs WHERE status='queued' ORDER BY created LIMIT 1").fetchone()
                if row:
                    con.execute("UPDATE jobs SET status='running',updated=? WHERE id=?",(time.time(),row['id']))
                    (root()/row['id']/'.worker-running').touch()
            if not row:await asyncio.sleep(1);continue
            jid=row['id'];folder=root()/jid
            try:
                with (folder/'worker.log').open('ab') as log:
                    process=await asyncio.create_subprocess_exec(sys.executable,'-m','paper2smiles.worker',jid,cwd=ROOT,stdout=log,stderr=log,start_new_session=True)
                    deadline=time.monotonic()+12*3600
                    while process.returncode is None:
                        try:await asyncio.wait_for(process.wait(),timeout=1)
                        except asyncio.TimeoutError:pass
                        with connect() as con:state=con.execute('SELECT status FROM jobs WHERE id=?',(jid,)).fetchone()[0]
                        if state=='cancelled' or time.monotonic()>deadline:
                            if process.returncode is None:os.killpg(process.pid,signal.SIGTERM)
                            try:await asyncio.wait_for(process.wait(),10)
                            except asyncio.TimeoutError:os.killpg(process.pid,signal.SIGKILL);await process.wait()
                            break
                with connect() as con:
                    state=con.execute('SELECT status FROM jobs WHERE id=?',(jid,)).fetchone()[0]
                    if state!='cancelled':
                        completed=process.returncode==0 and (folder/'report/report-data.json').exists()
                        message='提取未完成，可重试并继续检查点。'
                        try:
                            failure=json.loads((folder/'status.json').read_text()).get('error','')
                            if 'ACCOUNT_BALANCE_ERROR' in failure:message='模型 API 余额不足，请联系管理员充值后重试。'
                            elif 'model_not_found' in failure:message='模型 API 暂无可用通道，请联系管理员。'
                        except (OSError,ValueError):pass
                        con.execute('UPDATE jobs SET status=?,updated=?,error=? WHERE id=?',('completed' if completed else 'failed',time.time(),'' if completed else message,jid))
            except Exception:
                with connect() as con:con.execute("UPDATE jobs SET status='failed',updated=?,error=? WHERE id=? AND status!='cancelled'",(time.time(),'任务执行失败，请重试。',jid))
            finally:
                if process and process.returncode is None:
                    os.killpg(process.pid,signal.SIGTERM)
                    try:await asyncio.wait_for(process.wait(),10)
                    except asyncio.TimeoutError:os.killpg(process.pid,signal.SIGKILL);await process.wait()
                process=None
                (folder/'.worker-running').unlink(missing_ok=True)
    finally:
        if process and process.returncode is None:
            os.killpg(process.pid,signal.SIGTERM)
            try:await asyncio.wait_for(process.wait(),10)
            except asyncio.TimeoutError:os.killpg(process.pid,signal.SIGKILL);await process.wait()

@asynccontextmanager
async def lifespan(app):
    task=asyncio.create_task(worker_loop())
    yield
    task.cancel()
    try:await task
    except asyncio.CancelledError:pass

async def error(request,exc):return JSONResponse({'error':exc.detail},status_code=exc.status_code)
app=Starlette(routes=[Route('/health',health),Route('/jobs',jobs,methods=['GET','POST']),Route('/jobs/{job_id}',detail),Route('/jobs/{job_id}/cancel',cancel,methods=['POST']),Route('/jobs/{job_id}/retry',retry,methods=['POST']),Route('/jobs/{job_id}/files/{path:path}',artifact)],lifespan=lifespan,exception_handlers={HTTPException:error})
