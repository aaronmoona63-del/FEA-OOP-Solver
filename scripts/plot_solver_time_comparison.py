import matplotlib.pyplot as plt
import numpy as np
import os

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman', 'DejaVu Serif']
plt.rcParams['axes.unicode_minus'] = False    
plt.rcParams['font.size'] = 12

time_costs = [0.0, 0.0, 0.0]
if os.path.exists('benchmark_time.txt'):
    with open('benchmark_time.txt', 'r') as f:
        lines = f.readlines()
        time_costs[0] = float(lines[0].strip())
        time_costs[1] = float(lines[1].strip())
        time_costs[2] = float(lines[2].strip())
else:
    print("❌ 找不到 benchmark_time.txt！")
    exit()

solvers = ['Intel MKL Pardiso\n(Direct Solver)', 'Standard CG\n(No Preconditioner)', 'Jacobi PCG\n(Proposed Method)']
# 注意：这里的迭代次数是示意，写论文时请根据终端里显示的真实迭代次数在这行修改
iterations = [1, 722, 435] 

fig, ax1 = plt.subplots(figsize=(10, 6), dpi=150)
colors = ['#d62728', '#7f7f7f', '#1f77b4'] 
bars = ax1.bar(solvers, time_costs, color=colors, edgecolor='black', linewidth=1.5, width=0.5)

ax1.set_title('Computational Time Comparison of Different Linear Solvers\n(Regular Hexahedral Mesh, ~31k DOF)', fontsize=16, fontweight='bold', pad=20)
ax1.set_ylabel('Wall-clock Time (Seconds)', fontsize=14, fontweight='bold')

# 🌟 修复点：去掉了这里会报错的 fontweight='bold'
ax1.tick_params(axis='y', labelsize=12)
ax1.tick_params(axis='x', labelsize=13) 
ax1.grid(axis='y', linestyle='--', alpha=0.6)

for bar, iters in zip(bars, iterations):
    height = bar.get_height()
    ax1.text(bar.get_x() + bar.get_width()/2., height + (max(time_costs)*0.02), f'{height:.3f} s', ha='center', va='bottom', fontsize=13, fontweight='bold')
    iter_text = "Direct\n(No Iter)" if (iters == 1 or iters == 0) else f"{iters} Iters"
    ax1.text(bar.get_x() + bar.get_width()/2., height / 2, iter_text, ha='center', va='center', fontsize=12, color='white', fontweight='bold', bbox=dict(facecolor='black', alpha=0.5, edgecolor='none', pad=2))

ax1.set_ylim(0, max(time_costs) * 1.25)
plt.tight_layout()
plt.savefig('solver_time_comparison.png', dpi=300)
print("✅ 柱状图已成功生成并保存为 solver_time_comparison.png")
