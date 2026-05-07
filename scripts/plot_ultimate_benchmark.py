import matplotlib.pyplot as plt
import numpy as np
import os

# 学术字体安全设置
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman', 'DejaVu Serif']
plt.rcParams['axes.unicode_minus'] = False    

# 自动读取规则板耗时
time_regular = [0.0, 0.0, 0.0]
if os.path.exists('benchmark_time.txt'):
    with open('benchmark_time.txt', 'r') as f:
        lines = f.readlines()
        time_regular = [float(lines[0].strip()), float(lines[1].strip()), float(lines[2].strip())]

# 自动读取带孔板耗时
time_hole = [0.0, 0.0, 0.0]
if os.path.exists('benchmark_time_hole.txt'):
    with open('benchmark_time_hole.txt', 'r') as f:
        lines = f.readlines()
        time_hole = [float(lines[0].strip()), float(lines[1].strip()), float(lines[2].strip())]

# ==============================================================
# 🚨 唯一的填空区：请把刚才控制台打印出来的真实迭代次数填在这里 🚨
# ==============================================================
iter_regular = [0, 722, 435]       # 填入规则板的: [0, CG迭代次数, PCG迭代次数]
iter_hole = [0, 1614, 841]        # 填入带孔板的: [0, CG迭代次数, PCG迭代次数]
# ==============================================================

labels = ['Intel MKL Pardiso\n(Direct Solver)', 'Standard CG\n(No Preconditioner)', 'Jacobi PCG\n(Proposed Method)']
x = np.arange(len(labels))
width = 0.35  

fig, ax = plt.subplots(figsize=(11, 7), dpi=150)

rects1 = ax.bar(x - width/2, time_regular, width, label='Regular Mesh (Well-conditioned)', color='#1f77b4', edgecolor='black', linewidth=1.2)
rects2 = ax.bar(x + width/2, time_hole, width, label='Distorted Mesh with Hole (Ill-conditioned)', color='#d62728', edgecolor='black', linewidth=1.2)

ax.set_ylabel('Wall-clock Time (Seconds)', fontsize=14)
ax.set_title('Robustness Benchmark: Regular vs. Distorted Mesh (~31k DOF)', fontsize=16, pad=20)
ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=12)
ax.legend(fontsize=12, loc='upper left')
ax.grid(axis='y', linestyle='--', alpha=0.6)
ax.tick_params(axis='y', labelsize=12)
ax.tick_params(axis='x', labelsize=12)

def autolabel(rects, iters, times_arr):
    for rect, it in zip(rects, iters):
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., height + 0.02 * max(max(time_regular), max(time_hole)),
                f'{height:.3f} s', ha='center', va='bottom', fontsize=11)
        text = "Direct\n(No Iter)" if it == 0 else f"{it} Iters"
        ax.text(rect.get_x() + rect.get_width()/2., height / 2, text, 
                ha='center', va='center', fontsize=11, color='white',
                bbox=dict(facecolor='black', alpha=0.5, edgecolor='none', pad=2))

autolabel(rects1, iter_regular, time_regular)
autolabel(rects2, iter_hole, time_hole)

ax.set_ylim(0, max(max(time_regular), max(time_hole)) * 1.3)
plt.tight_layout()
plt.savefig('ultimate_robustness_benchmark.png', dpi=300)
print("✅ 终极双柱状图已成功生成并保存为 ultimate_robustness_benchmark.png")
