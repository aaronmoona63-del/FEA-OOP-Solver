import matplotlib.pyplot as plt
import numpy as np

# 设置纯净学术字体
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman', 'DejaVu Serif']
plt.rcParams['axes.unicode_minus'] = False    

# ==============================================================
# 🎯 严谨的 600k 规模实测数据 (解决标签重叠版)
# ==============================================================
dof = [31800, 222480, 422400, 588120]
time_pardiso = [0.612, 38.416, 140.658, 286.910]   
time_cg      = [3.100, 60.837, 142.119, 158.994]  
time_pcg     = [2.170, 86.623, 190.299, 297.385]  

fig, ax = plt.subplots(figsize=(11, 7.5), dpi=150)

# 绘制曲线
ax.plot(dof, time_pardiso, marker='o', markersize=9, linewidth=2.5, color='#E15759', label='Intel MKL Pardiso (Direct Method)')
ax.plot(dof, time_cg, marker='s', markersize=9, linewidth=3.0, color='#4E79A7', label='Standard CG (Iterative)')
ax.plot(dof, time_pcg, marker='^', markersize=9, linewidth=2.5, linestyle='--', color='#76B7B2', label='Jacobi PCG (Preconditioned)')

# 为每一个标签定制绝对精确的 (X, Y) 偏移量，彻底告别重叠
# 格式: (X偏移, Y偏移)
pardiso_offsets = [(0, -18), (0, -18), (15, -18), (-20, -15)]
pardiso_ha      = ['center', 'center', 'left', 'right']

cg_offsets      = [(0, 15),  (25, -5), (-15, 15), (0, -22)]
cg_ha           = ['center', 'left', 'right', 'center']

pcg_offsets     = [(25, -5), (0, 15),  (0, 15),   (15, 15)]
pcg_ha          = ['left', 'center', 'center', 'left']

# 循环写入时间标签
for i in range(4):
    ax.annotate(f"{time_pardiso[i]:.1f}s", (dof[i], time_pardiso[i]), 
                textcoords="offset points", xytext=pardiso_offsets[i], ha=pardiso_ha[i], 
                fontsize=11, color='#E15759', fontweight='bold')
    
    ax.annotate(f"{time_cg[i]:.1f}s", (dof[i], time_cg[i]), 
                textcoords="offset points", xytext=cg_offsets[i], ha=cg_ha[i], 
                fontsize=11, color='#4E79A7', fontweight='bold')
                
    ax.annotate(f"{time_pcg[i]:.1f}s", (dof[i], time_pcg[i]), 
                textcoords="offset points", xytext=pcg_offsets[i], ha=pcg_ha[i], 
                fontsize=11, color='#76B7B2', fontweight='bold')

# 核心标注 1：性能交叉区 (CG 打败 Pardiso)
ax.axvspan(380000, 450000, color='#F2CF5B', alpha=0.2)
ax.annotate('Scalability Crossover\n(Iterative overtakes Direct)', xy=(422400, 141), xytext=(200000, 220),
            arrowprops=dict(facecolor='#B8860B', shrink=0.05, width=1.5, headwidth=6),
            fontsize=12, color='#B8860B', fontweight='bold')

# 核心标注 2：解释 PCG 为何变慢 
ax.annotate('Hardware Overhead > Iteration Gain\n(Preconditioner matrix operations dominate)', 
            xy=(588120, 297.385), xytext=(350000, 320),
            arrowprops=dict(facecolor='#76B7B2', shrink=0.05, width=1.5, headwidth=6),
            fontsize=11, color='#76B7B2', bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#76B7B2", alpha=0.8))

# 设置坐标轴与外观
ax.set_xlabel('Total Degrees of Freedom (DOF)', fontsize=14, fontweight='bold')
ax.set_ylabel('Wall-clock Time (Seconds)', fontsize=14, fontweight='bold')
ax.set_title('Solver Scalability Comparison (Up to 600k DOF)', fontsize=16, fontweight='bold', pad=20)

xticks = np.arange(0, 700000, 100000)
ax.set_xticks(xticks)
ax.set_xticklabels([f"{int(x/1000)}k" for x in xticks], fontsize=12)
ax.set_xlim(0, 650000)
ax.set_ylim(0, 360)
ax.tick_params(axis='both', labelsize=12)
ax.legend(fontsize=11, loc='upper left')
ax.grid(True, linestyle='--', alpha=0.5)

plt.tight_layout()
plt.savefig('scalability_curve_600k_noclash.png', dpi=300)
print("✅ 排版完美的无重叠实测图表已生成！")
