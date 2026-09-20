# @ModuleName: Config
# @Function: 
# @Author: ggl
# @Time: 2026/7/21 17:00
# ==========================================================
# NCBI Configuration
# ==========================================================

# 建议填写自己的邮箱（NCBI官方要求）
EMAIL = "your_email@example.com"
# 如果以后申请了NCBI API Key，可以填这里（没有就留空）
NCBI_API_KEY = ""

# 每次请求之间暂停（秒）
# 没有API Key建议0.35左右
REQUEST_DELAY = 0.35

# 查询失败自动重试次数
MAX_RETRY = 3

# 输出编码
OUTPUT_ENCODING = "utf-8-sig"

# SQLite缓存数据库
CACHE_DB = "biosample_cache.db"

# 日志文件
LOG_FILE = "query.log"

# 失败accession保存
FAILED_FILE = "failed_accession.txt"

# 汇总统计
SUMMARY_FILE = "summary.xlsx"



# 批量 AI 处理时每批最多样本数（防止 token 超限）
AI_BATCH_SIZE = 50
# ==========================================================
# AI Configuration
# ==========================================================

# 是否启用AI分类
USE_AI = True

# OpenAI兼容接口
API_BASE = "https://api.deepseek.com"

# API Key
DEEPSEEK_API_KEY  = "sk-bed9399fc12b4d81abf2b0ad54356a02"

# 模型
MODEL = "deepseek-chat"