# 0.1.1

- 自动输出 00_results_summary.tsv 与 00_results_report.html。
- 新增 summarize --run：从现有统一版本结果生成简报，不重跑 BLAST。
- 弱或片段候选不进入主要系统候选列；同一 anchor 的竞争类型用斜线显示。
- 失败样本单独标记，未检出限定于当前库和阈值。
- 保护数据库原始字节及 SHA256 清单，增加独立 wheel 安装校验脚本。
- 提供 Conda 配方、构建和显式上传脚本，声明 BLAST 运行依赖。
- 检测阈值、数据库内容及功能证据等级未改变。

Conda 构建与渠道上传尚需在有 Conda 和发布账号的环境执行。
