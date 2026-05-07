import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from PIL import Image, ImageChops
import os

# ==========================================
# 1. 设置中文字体 (Windows 专属，防止中文显示为方块)
# ==========================================
plt.rcParams['font.sans-serif'] = ['SimHei']  
plt.rcParams['axes.unicode_minus'] = False 

# ==========================================
# 2. 自动裁剪白边函数 (保留核心黑科技)
# ==========================================
def autocrop_image(image_path):
    im = Image.open(image_path).convert("RGB")
    bg = Image.new(im.mode, im.size, (255, 255, 255))
    diff = ImageChops.difference(im, bg)
    diff = ImageChops.add(diff, diff, 2.0, -100)
    bbox = diff.getbbox()
    if bbox:
        pad = 20
        bbox = (max(0, bbox[0]-pad), max(0, bbox[1]-pad), 
                min(im.width, bbox[2]+pad), min(im.height, bbox[3]+pad))
        return im.crop(bbox)
    return im

# ==========================================
# 3. 图片配置与排版基底
# ==========================================
images = ["4.1.1.1.png", "4.1.1.2.png", "4.1.1.3.png", "4.1.1.4.png", "4.1.1.5.png"]

labels = [
    "(a) 极粗网格 (40)", 
    "(b) 较粗网格 (500)", 
    "(c) 中等网格 (5000)", 
    "(d) 较细网格 (40000)", 
    "(e) 极细网格 (135000)"
]

fig = plt.figure(figsize=(16, 10))
# 因为去掉了坐标轴，减小了 hspace 使得上下两排更紧凑
gs = GridSpec(2, 6, figure=fig, wspace=0.1, hspace=0.2)

ax1 = fig.add_subplot(gs[0, 0:2])
ax2 = fig.add_subplot(gs[0, 2:4])
ax3 = fig.add_subplot(gs[0, 4:6])
ax4 = fig.add_subplot(gs[1, 1:3])
ax5 = fig.add_subplot(gs[1, 3:5])

axes = [ax1, ax2, ax3, ax4, ax5]

# ==========================================
# 4. 组装图像与文字 (纯净模式)
# ==========================================
for i, ax in enumerate(axes):
    if os.path.exists(images[i]):
        img = autocrop_image(images[i])
        ax.imshow(img)
    else:
        print(f"找不到图片 {images[i]}！")
    
    ax.axis('off') # 隐藏边框和自带坐标轴
    
    # 添加纯净的文字标注
    ax.set_title(labels[i], y=-0.12, fontsize=18, fontweight='bold')

plt.tight_layout()

# 导出超高清图
output_filename = "Figure_4-2_Mesh_Densities_Pure.png"
plt.savefig(output_filename, dpi=400, bbox_inches='tight', facecolor='white')
print(f"✅ 极致纯净的论文排版图已生成：{output_filename}")