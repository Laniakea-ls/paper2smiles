"""Run one extraction; API credentials stay in the private service configuration."""
import re,sys
from pathlib import Path
from paper2smiles.settings import ROOT,settings
sys.path.insert(0,str(ROOT/'tools/paper_si_workflow'))
from model import Luna
from workflow_v3 import ActivityWorkflow

def main():
    jid=sys.argv[1]
    if not re.fullmatch('[a-f0-9]{32}',jid):raise ValueError('Invalid job id')
    cfg=settings();folder=Path(cfg['jobs_root'])/jid;key=cfg['api_key']
    if not key:raise RuntimeError('Model API key not configured')
    flow=ActivityWorkflow(folder/'upload.pdf',folder,Luna(key,cfg['api_base'],cfg['model'],folder))
    try:flow.run()
    except Exception as exc:
        flow.status('提取未完成',failed=True,error=str(exc).replace(key,'[REDACTED]')[:1200]);return 1
    return 0
if __name__=='__main__':raise SystemExit(main())
