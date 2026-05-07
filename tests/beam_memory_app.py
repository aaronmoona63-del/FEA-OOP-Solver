import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# ==========================================
# 页面基本设置
# ==========================================
st.set_page_config(page_title="3D FEA 内存与约束推演仪", layout="wide")
st.title("🧊 三维长方形梁：网格、约束与稀疏矩阵内存推演仪")
st.markdown("通过可视化三维实体网格与边界约束，深度剖析不同矩阵存储格式的空间复杂度。")

# ==========================================
# 侧边栏：网格与约束控制
# ==========================================
st.sidebar.header("⚙️ 物理模型参数设置")
st.sidebar.markdown("模拟 Hex8 八节点六面体单元构成的悬臂梁/简支梁")

# 滑块设置 (为了保证 3D 绘图不卡顿，此处规模限制在教学级，但图表推演到工业级)
nx = st.sidebar.slider("X方向单元数 (Nx)", min_value=1, max_value=20, value=5, step=1)
ny = st.sidebar.slider("Y方向单元数 (Ny)", min_value=1, max_value=20, value=5, step=1)
nz = st.sidebar.slider("Z方向单元数 (Nz, 梁长)", min_value=2, max_value=50, value=20, step=1)

bc_type = st.sidebar.selectbox(
    "🔒 边界约束条件 (Dirichlet BC)",
    ("无约束 (Free)", "底面全固支悬臂梁 (Z=0 锁定)", "两端全固支 (Z=0 与 Z=Nz 锁定)")
)

# ==========================================
# 核心力学与数学逻辑分析
# ==========================================
# 1. 节点计算逻辑
nodes_x, nodes_y, nodes_z = nx + 1, ny + 1, nz + 1
total_nodes = nodes_x * nodes_y * nodes_z

# 2. 约束节点计算逻辑
constrained_nodes = 0
if bc_type == "底面全固支悬臂梁 (Z=0 锁定)":
    constrained_nodes = nodes_x * nodes_y
elif bc_type == "两端全固支 (Z=0 与 Z=Nz 锁定)":
    constrained_nodes = 2 * nodes_x * nodes_y

free_nodes = total_nodes - constrained_nodes
dof_N = free_nodes * 3  # 有效自由度 (方程组真实阶数)

# 3. 拓扑特征计算逻辑 (极其关键)
# 最大半带宽 B: 编号顺序为 X->Y->Z，最大跨度出现在 Z 方向相邻节点之间
bandwidth_B = 3 * (nodes_x * nodes_y) + 3

# 非零元素 NNZ 估算: 内部节点最多连接 27 个节点，每个连接是 3x3 块
nnz = dof_N * 81 if dof_N > 0 else 0

# Skyline 轮廓线长度积分近似
if dof_N > bandwidth_B:
    L_sky = dof_N * bandwidth_B - (bandwidth_B**2) / 2
else:
    L_sky = (dof_N**2) / 2

sparsity = (nnz / (dof_N**2)) * 100 if dof_N > 0 else 0

# ==========================================
# 模块 A：三维网格与约束可视化 (直观展示物理模型)
# ==========================================
st.subheader("📍 物理模型：长方体梁网格与约束状态图")

# 生成 3D 坐标点
X, Y, Z = np.meshgrid(np.linspace(0, nx, nodes_x), 
                      np.linspace(0, ny, nodes_y), 
                      np.linspace(0, nz, nodes_z), indexing='ij')

X_flat, Y_flat, Z_flat = X.flatten(), Y.flatten(), Z.flatten()

# 区分自由节点和约束节点
if bc_type == "底面全固支悬臂梁 (Z=0 锁定)":
    mask_fixed = (Z_flat == 0)
elif bc_type == "两端全固支 (Z=0 与 Z=Nz 锁定)":
    mask_fixed = (Z_flat == 0) | (Z_flat == nz)
else:
    mask_fixed = np.zeros_like(Z_flat, dtype=bool)

mask_free = ~mask_fixed

fig_3d = plt.figure(figsize=(10, 6))
ax_3d = fig_3d.add_subplot(111, projection='3d')

# 绘制散点
# 为了防止点太多浏览器卡死，做一个简单的采样显示
step = 1 if total_nodes < 8000 else 2 
ax_3d.scatter(X_flat[mask_free][::step], Y_flat[mask_free][::step], Z_flat[mask_free][::step], 
              c='#1f77b4', s=10, alpha=0.5, label='Free Nodes (Unconstrained)')
if np.any(mask_fixed):
    ax_3d.scatter(X_flat[mask_fixed], Y_flat[mask_fixed], Z_flat[mask_fixed], 
                  c='red', s=40, marker='s', label='Fixed Nodes (Dirichlet BC)')

ax_3d.set_xlabel('X')
ax_3d.set_ylabel('Y')
ax_3d.set_zlabel('Z (Length)')
ax_3d.set_title('3D Hex8 Beam Mesh & Boundary Conditions', fontweight='bold')
ax_3d.legend()
st.pyplot(fig_3d)

# ==========================================
# 模块 B：数据仪表盘与表格
# ==========================================
col1, col2, col3, col4 = st.columns(4)
col1.metric("网格总节点数", f"{total_nodes:,}")
col2.metric("被约束节点数", f"{constrained_nodes:,}")
col3.metric("有效总自由度 (DOF)", f"{dof_N:,}")
col4.metric("整体刚度矩阵稀疏度", f"{sparsity:.4f} %")

MB = 1048576.0
mem_dense = (dof_N**2 * 8) / MB
mem_banded = (dof_N * (2 * bandwidth_B + 1) * 8) / MB
mem_skyline = (L_sky * 8 + (dof_N + 1) * 4) / MB
mem_coo = (nnz * 16) / MB
mem_csr = (nnz * 12 + (dof_N + 1) * 4) / MB

st.markdown("### 📊 当前物理状态下内存消耗对比 (MB)")
data = {
    "存储格式": ["Dense (全矩阵)", "Banded (等带宽)", "Skyline (轮廓线)", "COO (坐标系)", "CSR (压缩稀疏行)"],
    "理论空间复杂度": ["O(N²)", "O(N * B)", "O(N * B) 优化", "O(NNZ)", "O(NNZ) 极限压缩"],
    "消耗内存 (MB)": [mem_dense, mem_banded, mem_skyline, mem_coo, mem_csr]
}
df = pd.DataFrame(data)
st.dataframe(df.style.format({"消耗内存 (MB)": "{:,.2f}"}), use_container_width=True)

# ==========================================
# 模块 C：理论内存推演图 (Log-Log)
# ==========================================
st.markdown("### 📈 理论内存消耗随自由度演化规律 (用于论文验证)")

N_arr = np.logspace(3, 6, 100) # 从 1000 到 100万 自由度
B_arr = 3 * (N_arr/3)**(2/3) + 3 # 近似 3D 带宽
NNZ_arr = N_arr * 81
L_sky_arr = N_arr * B_arr - (B_arr**2)/2

fig_log, ax_log = plt.subplots(figsize=(10, 5))
ax_log.loglog(N_arr, (N_arr**2 * 8)/MB, 'r-', linewidth=2, label='Dense')
ax_log.loglog(N_arr, (N_arr * (2*B_arr + 1) * 8)/MB, 'orange', linestyle='--', linewidth=2, label='Banded')
ax_log.loglog(N_arr, (L_sky_arr * 8 + (N_arr + 1) * 4)/MB, 'm-.', linewidth=2, label='Skyline')
ax_log.loglog(N_arr, (NNZ_arr * 16)/MB, 'g:', linewidth=2, label='COO')
ax_log.loglog(N_arr, (NNZ_arr * 12 + (N_arr + 1) * 4)/MB, 'b-', linewidth=3, label='CSR (OOFEM Solver)')

# 标出当前模型所处的点
ax_log.scatter([dof_N]*5, [mem_dense, mem_banded, mem_skyline, mem_coo, mem_csr], 
               color='black', zorder=5, s=60, edgecolors='white', label='Current Model State')

ax_log.set_xlabel('Effective Degrees of Freedom (N)', fontweight='bold')
ax_log.set_ylabel('Memory Consumption (MB)', fontweight='bold')
ax_log.set_title('Asymptotic Space Complexity Comparison', fontweight='bold')
ax_log.grid(True, which="both", ls="--", alpha=0.5)
ax_log.legend()
st.pyplot(fig_log)