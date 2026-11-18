# Paper2SMILES

文献活性提取后端。上传 PDF 后，先确定具有数值或测定限值的活性编号，再按作者 SI → 原文完整名称 / OPSIN → 通式连接 → 完整结构图三模型识别的顺序寻找结构。完整结构图的三模型预测中，两路规范化 SMILES 一致即可采用；无共识的结果保留待核对。

本仓库只维护 Python API、提取流程、报告数据、结构图片与 Excel 导出。前端在 [aidd_portal_frontend](https://github.com/Laniakea-ls/aidd_portal_frontend)：`src/tools/paper-sar` 复用 SAR 的 WorkflowShell、上传控件和结果组件；`src/tools/literature-report` 渲染报告。

## 运行

需要 Python 3.11+，以及本机已部署的 SAR 裁图 / PaddleOCR 服务、MolParser-Mobile-V2、MolSight-stereo、OCSRGlyph 和 OPSIN。现有模型适配器默认使用 AIDD 服务器目录，迁移服务器时需调整这些适配路径。仓库不包含模型权重、用户 PDF、识别缓存或密钥。

```sh
python3.11 -m venv .venv
.venv/bin/pip install -e .
cp config/service.example.json config/service.local.json
chmod 600 config/service.local.json
# 填写模型 API、至少 32 字符随机 identity_secret 和任务目录
.venv/bin/python -m uvicorn paper2smiles.server:app --host 127.0.0.1 --port 8118 --workers 1
```

`PAPER2SMILES_HOME` 指定运行数据根目录（默认当前工作目录）；`PAPER2SMILES_CONFIG` 可指定配置路径；`PAPER_SI_API_KEY` 可覆盖模型密钥。模型默认 `GPT-6-Luna`，标识区分大小写。部署可使用 `deploy/paper2smiles.service`，先调整工作目录与 Python 路径。

Portal 使用现有 `tool_gateway` 注册 `config/portal-tool.example.json` 中的工具配置。网关验证登录后签发绑定用户、请求方法和路径的短期 HS256 JWT；后端仅监听回环地址。前端通过 `/api/paper-sar` 访问后端，不持有模型或服务密钥。

## 接口

- `GET /health`：配置状态。
- `POST /jobs?filename=article.pdf`：以 `application/pdf` 上传 PDF，返回任务编号；最多 50 MB / 500 页。
- `GET /history`：门户统一历史记录（按登录用户隔离，包含状态、进度、报告及 Excel 操作）。
- `GET /jobs`、`GET /jobs/{id}`：当前用户的任务列表、阶段、汇总和报告路径。
- `POST /jobs/{id}/cancel`、`POST /jobs/{id}/retry`：取消或从已有检查点重试。
- `GET /jobs/{id}/files/report-data.json`、`assets/*`、`downloads/*`：已完成任务的报告、图片和 Excel。

SQLite 保存队列，单 worker 控制模型资源占用，重启后恢复未完成任务。所有接口均要求网关签名，任务及文件按用户隔离。相同用户重复上传同一 PDF 返回已有运行中 / 已完成任务。

Excel 以文章标题命名，包含化合物编号、SMILES、2D 结构、实验条件、合并的活性列和出处。正文与 SI 原始数据保留来源区别。无法确认的结构不补造；三模型一致仅代表识别共识，不代表独立准确率。

## 仓库结构

```text
src/paper2smiles/    Python 包
  server.py         应用工厂与依赖组装
  api.py            ExtractionAPI：请求处理
  authentication.py GatewayAuthenticator：网关身份验证
  jobs/             JobRepository、JobPresenter、JobRunner
  worker.py         ExtractionWorker：单任务执行
  extraction/       提取流程、结构解析与 ReportBuilder
tests/             unit/ 提取测试；service/ API 与队列测试
scripts/experiments/ 历史实验脚本，不随包安装
config/             配置示例（本地配置不入库）
deploy/             systemd 部署配置
```

流程由 `ActivityWorkflow` 管理，报告由 `ReportBuilder` 生成；服务通过应用工厂注入配置，每个实例拥有独立依赖。无状态的化学解析和数值转换保留为函数。Python 使用 PEP 8 风格、79 列格式与 NumPy 风格文档字符串。

## 开发与检查

```sh
.venv/bin/pip install -e '.[dev]'
.venv/bin/ruff check .
.venv/bin/ruff format --check .
.venv/bin/python -m unittest discover -s tests -v
.venv/bin/python -m build
```

测试需要上述部署环境中的化学依赖。`src` 布局要求先安装包；请勿通过修改 `sys.path` 导入源码。运行数据 `runs/`、本地配置、模型权重和构建产物不提交。

安装后也可使用 `paper2smiles-extract`、`paper2smiles-batch` 和 `paper2smiles-table-ocr` 命令；使用 `--help` 查看参数。
