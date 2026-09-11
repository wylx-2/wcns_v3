# v1.1.0 发布后外部 CI 移除记录

状态：**外部 CI 已移除，本机回归和 clean-HEAD 源码包卡口全部通过。**

## 1. 决策依据

v1.1.0 首次向 GitHub 同步 `main`、release 和多个阶段分支时，`wcns-ci` 的 `push` 触发器
同时启动 9 次运行。[示例运行 `34615909522`](https://github.com/wylx-2/wcns_v3/actions/runs/34615909522)
包含 5 个作业：Linux GCC、Linux Clang 和 Clang ASan+UBSan 成功，Linux OpenMPI 测试步骤和
Windows MSVC 配置步骤失败。这些外部结果
不属于已批准的本机版本卡口，持续触发只会产生无效告警和资源消耗。

项目负责人确认预期开发继续使用本机独立流程，不使用外部 CI。因此删除 GitHub Actions
工作流比继续修复或静默忽略外部作业更符合当前维护边界。

## 2. 移除范围

- 删除唯一的 `.github/workflows/ci.yml`，停止新分支从当前 `main` 继承 `wcns-ci`。
- 从 `tools/package_release.py` 删除对该工作流文件的强制打包依赖。
- 更新 README、运行/开发指南、已知限制、路线图、版本计划及发布说明中的当前状态。
- 保留 P/U 阶段设计与验收文档对历史工作流的事实记录，并增加发布后说明，不改写既有证据。

## 3. 明确保留

本次不删除或弱化以下本机质量能力：

- CMake 的 `WCNS_BUILD_TESTS`、`WCNS_ENABLE_MPI`、`WCNS_ENABLE_CGNS` 等构建选项；
- 串行/MPI CTest、算法规格检查、发布矩阵、失败路径和算例验证脚本；
- clang-format、Black、Clang-Tidy/Pylint 等本机格式与静态检查方法；
- v1.1.0 正式标签、历史 Actions 运行及阶段验收记录。

## 4. 自动验收

| 卡口 | 结果 |
|---|---:|
| 仓库工作流与当前 CI 引用审计 | 通过；唯一工作流已删除，历史引用已标注 |
| Python 格式、编译与打包器定向检查 | 32/32 Black 通过；打包器 `py_compile` 通过 |
| 算法规格检查 | 6/6 通过 |
| 本地 CTest 回归 | 串行 60/60（26.05 s），MPI 108/108（77.55 s） |
| 源码包双次生成、逐字节比较与内容审计 | 通过；341 个文件，无 workflow，Case07 保留，NACA0012 排除 |
| Markdown 相对链接及 `git diff --check` | 65 份文档通过；diff 通过 |

正式 `v1.1.0` 标签保持不可移动。本维护提交位于标签之后，只改变开发基础设施和文档，不改变
求解器源码、数值算法、配置语义、测试内容或已有计算结果。
