# CGNS dependency

`CGNS-4.4.0.zip` 复制自 PHengLEI 随附的 CGNS 4.4.0 源码包：

```text
PHengLEI/PHengLEI/3rdparty/CGNS-4.4.0.zip
```

SHA-256：

```text
72AD499936089102C1D77E3CDD78A623D61D1C6AD94C4DF3B23F36DB9DDFAD6A
```

WCNS 当前使用静态 ADF 后端构建，关闭 HDF5、Fortran、共享库和 CGNS 自带测试。这样可以在当前 MinGW 工具链上构建，同时避免依赖 PHengLEI 中由 MSVC 生成、无法与 MinGW 稳定混用的 `cgns.lib`。

