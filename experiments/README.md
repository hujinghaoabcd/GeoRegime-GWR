# Experiments

本目录用于论文实验脚本。不要把一次性的 notebook 逻辑直接当成最终可复现实验。

建议子目录：

- `simulations/`：合成 DGP 与理论验证；
- `ablations/`：组件消融；
- `benchmarks/`：与 GWR / spatial regime / clustered varying coefficient 等方法对比；
- `case_studies/`：真实数据案例；
- `figures/`：论文图生成脚本。

每个正式实验至少记录 seed、参数、数据配置、模型配置、commit SHA 和原始指标。
