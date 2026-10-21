# Conda 构建与发布 0.1.1

目前仓库提供配方和脚本；不代表已经在 Anaconda.org 或 Bioconda 上架。
Python >=3.9，BLAST >=2.16,<3；无需 Prokka、ABRicate、Biopython。
软件为 noarch Python 包，但完整运行还要求目标平台有可用 BLAST 包；优先验证 Linux x86_64。

## 开发者构建（Linux）

```bash
conda create -n cbf-build --override-channels -c conda-forge python=3.11 conda-build anaconda-client pip setuptools wheel -y
conda activate cbf-build
cd /path/to/ColiBacFinder
bash scripts/check_package.sh
bash scripts/build_conda.sh
```

脚本不自动上传。conda build 会在测试环境验证 CLI、数据库校验和 BLAST 依赖。
产物在 conda-dist/noarch；若出错，停止发布并保留日志。

## 开发者上传

先自行注册 Anaconda.org 账号；该账号不一定与 GitHub 用户名相同。
不应把密码或 token 写进代码或发到聊天。

```bash
anaconda login
bash scripts/publish_conda.sh YOUR_ANACONDA_USERNAME /absolute/path/to/colibacfinder-0.1.1-py_0.conda
```

产物后缀和 build 字符串以实际生成文件为准（可能是 .tar.bz2）。
上传前确认源码及第三方参考数据的分发许可；本仓库尚未声明开源许可证。
脚本会重新测试包，上传不使用 --force。只有上传成功并确认渠道可见后才能分发下面的安装命令。

## 用户安装（发布完成后）

```bash
conda create -n colibacfinder --override-channels -c YOUR_ANACONDA_USERNAME -c conda-forge -c bioconda colibacfinder=0.1.1 -y
conda activate colibacfinder
colibacfinder db-check
colibacfinder run -i genomes -o new_results -t 8
```

这是个人 Conda 渠道发布，不是已被 Bioconda 接收。用户安装无需 pip 或开发构建工具。
