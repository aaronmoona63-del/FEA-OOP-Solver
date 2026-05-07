import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ==========================================
# 页面基本设置
# ==========================================
st.set_page_config(page_title="3D FEA 内存推演仪 V2", layout="wide")
st.title("🧊 三维长方形梁：高规模网格存储复杂度推演仪")
st.markdown("本工具支持高密度网格模拟，用于论证 **CSR** 格式在解决大规模三维实体问题时的核心价值。")

# ==========================================
# 侧边栏：网格与约束控制
# ==========================================
st.sidebar.header("⚙️ 物理模型参数设置")
st.sidebar.markdown("设定三维 Hex8 实体网格规模")

# 按照用户要求更新上限
nx = st.sidebar.slider("X方向单元数 (Nx)", min_value=1, max_value=30, value=10, step=1)
ny = st.sidebar.slider("Y方向单元数 (Ny)", min_value=1, max_value=30, value=10, step=1)
nz = st.sidebar.slider("Z方向单元数 (Nz, 梁长)", min_value=2, max_value=150, value=50, step=5)

bc_type = st.sidebar.selectbox(
    "🔒 边界约束条件 (Dirichlet BC)",
    ("无约束 (Free)", "底面全固支悬臂梁 (Z=0 锁定)", "两端全固支 (Z=0 与 Z=Nz 锁定)")
)

# ==========================================
# 核心力学与数学逻辑分析
# ==========================================
nodes_x, nodes_y, nodes_z = nx + 1, ny + 1, nz + 1
total_nodes = nodes_x * nodes_y * nodes_z

# 计算约束
constrained_nodes = 0
if bc_type == "底面全固支悬臂梁 (Z=0 锁定)":
    constrained_nodes = nodes_x * nodes_y
elif bc_type == "两端全固支 (Z=0 与 Z=Nz 锁定)":
    constrained_nodes = 2 * nodes_x * nodes_y

free_nodes = total_nodes - constrained_nodes
dof_N = free_nodes * 3  # 每个节点 3 个位移自由度

# 拓扑特征计算 (基于 3D 规则网格编号)
bandwidth_B = 3 * (nodes_x * nodes_y) + 3
nnz = dof_N * 81 if dof_N > 0 else 0 # 理论估算：每个节点及其邻居的 3x3 块

# Skyline 轮廓线长度
if dof_N > bandwidth_B:
    L_sky = dof_N * bandwidth_B - (bandwidth_B**2) / 2
else:
    L_sky = (dof_N**2) / 2

sparsity = (nnz / (dof_N**2)) * 100 if dof_N > 0 else 0

# ==========================================
# 模块 A：数据可视化与 3D 示意
# ==========================================
col_metrics, col_plot = st.columns([1, 2])

with col_metrics:
    st.metric("总自由度 (DOF)", f"{dof_N:,}")
    st.metric("非零元个数 (NNZ)", f"{nnz:,}")
    st.metric("矩阵稀疏度", f"{sparsity:.5f} %")
    
    # 内存计算 (MB)
    MB = 1048576.0
    mem_dense = (dof_N**2 * 8) / MB
    mem_banded = (dof_N * (2 * bandwidth_B + 1) * 8) / MB
    mem_skyline = (L_sky * 8 + (dof_N + 1) * 4) / MB
    mem_coo = (nnz * 16) / MB
    mem_csr = (nnz * 12 + (dof_N + 1) * 4) / MB

with col_plot:
    # 3D 渲染预览 (带采样保护)
    fig_3d = plt.figure(figsize=(8, 5))
    ax_3d = fig_3d.add_subplot(111, projection='3d')
    
    # 仅绘制边缘点以防卡顿
    if total_nodes < 5000:
        X, Y, Z = np.meshgrid(np.linspace(0, nx, nodes_x), np.linspace(0, ny, nodes_y), np.linspace(0, nz, nodes_z))
        ax_3d.scatter(X, Y, Z, c='blue', s=2, alpha=0.3)
    else:
        st.info("💡 当前网格规模巨大，3D 视图已自动切换为边界外廓显示模式以保证流畅。")
        # 简单绘制长方体线框
        for x_edge in [0, nx]:
            for y_edge in [0, ny]:
                ax_3d.plot([x_edge, x_edge], [y_edge, y_edge], [0, nz], 'gray', alpha=0.5)
        for z_edge in [0, nz]:
            for x_edge in [0, nx]:
                ax_3d.plot([x_edge, x_edge], [0, ny], [z_edge, z_edge], 'gray', alpha=0.5)
    
    ax_3d.set_title(f"Mesh: {nx}x{ny}x{nz} (Total Nodes: {total_nodes:,})")
    st.pyplot(fig_3d)

# ==========================================
# 模块 B：内存对比表格与 Log-Log 图
# ==========================================
st.markdown("### 📊 内存消耗对比与空间复杂度演化")

data = {
    "存储格式": ["Dense (全矩阵)", "Banded (等带宽)", "Skyline (轮廓线)", "COO (坐标系)", "CSR (压缩稀疏行)"],
    "内存消耗": [mem_dense, mem_banded, mem_skyline, mem_coo, mem_csr]
}
df = pd.DataFrame(data)
# 单位自动切换逻辑
def format_mem(x):
    if x > 1024 * 1024: return f"{x/(1024*1024):.2f} TB"
    if x > 1024: return f"{x/1024:.2f} GB"
    return f"{x:.2f} MB"

df["格式化显示"] = df["内存消耗"].apply(format_mem)
st.table(df[["存储格式", "格式化显示"]])

# 双对数图
N_arr = np.logspace(3, 6, 100)
B_arr = 3 * (N_arr/3)**(2/3) + 3
NNZ_arr = N_arr * 81
L_sky_arr = N_arr * B_arr

fig_log, ax_log = plt.subplots(figsize=(10, 5))
ax_log.loglog(N_arr, (N_arr**2 * 8)/MB, 'r-', label='Dense O(N²)')
ax_log.loglog(N_arr, (N_arr * (2*B_arr + 1) * 8)/MB, 'orange', ls='--', label='Banded O(N * N^(2/3))')
ax_log.loglog(N_arr, (L_sky_arr * 8)/MB, 'm-.', label='Skyline')
ax_log.loglog(N_arr, (NNZ_arr * 16)/MB, 'g:', label='COO O(NNZ)')
ax_log.loglog(N_arr, (NNZ_arr * 12)/MB, 'b-', lw=3, label='CSR (OOFEM Optimizer)')

ax_log.scatter([dof_N]*5, [mem_dense, mem_banded, mem_skyline, mem_coo, mem_csr], color='black', zorder=5)
ax_log.set_xlabel('DOF (N)')
ax_log.set_ylabel('Memory (MB)')
ax_log.grid(True, which="both", ls="--", alpha=0.5)
ax_log.legend()
st.pyplot(fig_log)