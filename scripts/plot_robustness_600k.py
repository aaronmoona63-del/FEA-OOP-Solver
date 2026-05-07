import matplotlib.pyplot as plt
import numpy as np

# 设置纯净的学术字体
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman', 'DejaVu Serif']
plt.rcParams['axes.unicode_minus'] = False    

# ==============================================================
# 🚨 你的全部真实数据 (约 60 万自由度级)
# ==============================================================
# 对照组：60万均匀板 (网格完美，条件数极佳)
time_regular = [175.29, 39.44, 53.82]   # [Pardiso, CG, PCG]
iter_regular = [0, 1291, 1002]

# 实验组：58.8万带孔板 (网格畸变，条件数恶劣)
time_hole = [286.91, 158.99, 297.38]    # [Pardiso, CG, PCG]
iter_hole = [0, 5193, 5527]
# ==============================================================

labels = ['Intel MKL Pardiso\n(Direct Solver)', 'Standard CG\n(No Preconditioner)', 'Jacobi PCG\n(Diagonal Precond.)']
x = np.arange(len(labels))
width = 0.35  

fig, ax = plt.subplots(figsize=(10.5, 7), dpi=150)

# 这次使用常规线性坐标轴即可，因为差距在同一数量级内，视觉冲击力最强
rects1 = ax.bar(x - width/2, time_regular, width, label='Uniform Plate (Ideal Condition Number)', color='#76B7B2', edgecolor='black', linewidth=1.2)
rects2 = ax.bar(x + width/2, time_hole, width, label='Plate with Hole (Poor Condition Number)', color='#E15759', edgecolor='black', linewidth=1.2)

ax.set_ylabel('Wall-clock Time (Seconds)', fontsize=14, fontweight='bold')
ax.set_title('Robustness Benchmark: Condition Number Impact on 600k DOF Systems', fontsize=16, fontweight='bold', pad=20)
ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=12, fontweight='bold')
ax.legend(fontsize=12, loc='upper left')
ax.grid(axis='y', linestyle='--', alpha=0.5)

# 自动标注时间和迭代次数
def autolabel(rects, times, iters):
    for rect, t, it in zip(rects, times, iters):
        height = rect.get_height()
        # 标注耗时 (秒)
        ax.text(rect.get_x() + rect.get_width()/2., height + 5,
                f'{t:.1f} s', ha='center', va='bottom', fontsize=12, fontweight='bold')
        
        # 标注迭代次数在柱子中间
        text = "Direct\n(No Iter)" if it == 0 else f"{it} Iters"
        ax.text(rect.get_x() + rect.get_width()/2., height / 2, text,
                ha='center', va='center', fontsize=11, color='white', fontweight='bold',
                bbox=dict(facecolor='black', alpha=0.6, edgecolor='none', boxstyle='round,pad=0.3'))

autolabel(rects1, time_regular, iter_regular)
autolabel(rects2, time_hole, iter_hole)

# 设置Y轴上限以留出文字空间
ax.set_ylim(0, max(max(time_regular), max(time_hole)) * 1.25)

plt.tight_layout()
plt.savefig('robustness_600k_benchmark.png', dpi=300)
print("✅ 终极 60万自由度鲁棒性对比大图已生成！请查看 robustness_600k_benchmark.png")
