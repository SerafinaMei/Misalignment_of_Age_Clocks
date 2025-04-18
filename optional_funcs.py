import os
import matplotlib.pyplot as plt

def plot_pie_chart(labels, sizes, colors, title, save_path=None):
    fig, ax = plt.subplots(figsize=(7, 5))
    wedges, texts, autotexts = ax.pie(
        sizes, labels=labels, autopct='%1.1f%%', colors=colors,
        startangle=140, wedgeprops={'edgecolor': 'white'}, textprops={'fontsize': 12}
    )
    ax.set_title(title, fontsize=14, fontweight='bold')
    
    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show()
