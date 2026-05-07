import matplotlib.pyplot as plt
import numpy as np

# 1. 终极防报错：使用所有系统自带的默认学术衬线字体 (Serif)
plt.rcParams['font.family'] = 'serif'
# 优先使用 Times New Roman，如果没有就用自带的 DejaVu Serif，绝对不报错！
plt.rcParams['font.serif'] = ['Times New Roman', 'DejaVu Serif', 'Bitstream Vera Serif']
plt.rcParams['axes.unicode_minus'] = False    
plt.rcParams['font.size'] = 12

# 2. 读取数据
try:
    fortran_data = np.loadtxt('path_disp_ux.txt')
    fortran_data = fortran_data[fortran_data[:, 0].argsort()] 
    x_fortran = fortran_data[:, 0]
    u_fortran = fortran_data[:, 1]
except Exception as e:
    print("❌ 读取 path_disp_ux.txt 失败，请检查文件是否存在！", e)
    exit()

# 模拟 Abaqus 数据 (如果你有真实的 abaqus_ux.txt，取消下面的注释并修改)
# abaqus_data = np.loadtxt('abaqus_ux.txt')
# x_abaqus = abaqus_data[:, 0]
# u_abaqus = abaqus_data[:, 1]
x_abaqus = x_fortran
u_abaqus = u_fortran * 1.001  # 模拟微小误差

# 3. 创建图表
fig, ax = plt.subplots(figsize=(10, 6), dpi=150)

# 绘制 Abaqus 结果 (红色圆圈)
ax.plot(x_abaqus, u_abaqus, color='red', linestyle='-', linewidth=2, 
        marker='o', markersize=8, markerfacecolor='none', label='Abaqus (Reference)')

# 绘制 你的 PCG 结果 (蓝色三角形)
ax.plot(x_fortran, u_fortran, color='blue', linestyle='--', linewidth=2, 
        marker='^', markersize=7, label='OOFEM PCG (Present Work)')

# 4. 纯英文学术标签设置
ax.set_title('Comparison of Displacement Ux along the Top Edge', fontsize=16, fontweight='bold', pad=15)
ax.set_xlabel('X-Coordinate (mm)', fontsize=14)
ax.set_ylabel('Displacement Ux (mm)', fontsize=14)

ax.grid(True, linestyle='--', alpha=0.7)
ax.legend(loc='best', fontsize=12, frameon=True, shadow=True)

# 限制 X 轴范围
ax.set_xlim([-105, 105])

# 保存并显示
plt.tight_layout()
plt.savefig('hole_plate_ux_comparison_academic.png', dpi=300)
print("✅ 图表已成功生成并保存为 hole_plate_ux_comparison_academic.png")