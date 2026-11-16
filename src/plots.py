import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

plt.rcParams.update({
    "font.size": 12,         # base size for all text
    "axes.labelsize": 14,   # Default: 14
    "xtick.labelsize": 12,  # Default: 12
    "ytick.labelsize": 12,  # Default: 12
    "legend.fontsize": 14,  # Default: 12
})

def sorted_xy(hw_to_count):
    """
        Hamming weights (x) and counts (y) of a histogram, sorted by Hamming
        weight.
    """
    xs = sorted(hw_to_count.keys())
    ys = [hw_to_count[x] for x in xs]
    return xs, ys


def save_pdf(fig, path, name):
    """
        Save fig as path/name.pdf.
    """
    fn = f'{path}/{name}.pdf'
    fig.savefig(fn, format="pdf", bbox_inches="tight", pad_inches=0)
    print(f"Wrote file {fn}")


def plot_dist(data_kernel, data_kernel_R, data_F3, expected_X_R, desc, path):
    """
        Generate plots for distribution of Z_i,
        (cf. Figure 2)
    """
    fig, ax = plt.subplots(figsize=(7, 2.8))

    # F3, kernel, and R lowest hw's of the kernel (averaged over experiments)
    x_F3, y_F3 = sorted_xy(data_F3)
    x_ker, y_ker = sorted_xy(data_kernel)
    x_R, y_R = sorted_xy(data_kernel_R)

    ax.bar(x_F3, y_F3, width=1, color='gray', alpha=.5,
           label='$\\mathbb{F}_{3}^{~l}$')
    ax.bar(x_ker, y_ker, width=.6, color='blue', alpha=.5,
           label='$\\text{hw}(\\vec{z}_1), ..., \\text{hw}(\\vec{z}_N)$')
    ax.bar(x_R, y_R, width=.6, color='purple', alpha=.8,
           label='$\\text{hw}(\\vec{z}_1), ..., \\text{hw}(\\vec{z}_R)$')

    ax.axvline(x=expected_X_R, color='red', linestyle='--', linewidth=1.2,
               label="$\\mathbb{E}[\\text{hw}(\\vec{z}_R)]$")

    # Setup
    xmin, xmax = ax.get_xlim()
    step = max(1, (xmax - xmin) // 5)
    ax.xaxis.set_major_locator(MultipleLocator(step))
    ax.set_ylabel("count")
    ax.legend(loc='upper right')

    save_pdf(fig, path, desc)

    # Same plot, log scale
    ax.set_yscale("log", base=2)
    ax.set_ylabel("log(count)")

    save_pdf(fig, path, desc + '_log')
    plt.close(fig)


def plot_attack(data_dump, desc, path):
    """
        Generate plots for demonstration of collision attack,
        (cf. Figure 3)
    """
    fig, ax = plt.subplots(figsize=(7, 2.8))

    success_prob = [y / data_dump['num_exps'] for y in data_dump['success']]

    ax.plot(data_dump['al'], success_prob, color='blue', marker='o',
            linestyle='-', alpha=.6, label='CollisionSearch')
    ax.axvline(x=data_dump['max_al'], color='red', linestyle='--',
               linewidth=1.2, label="prediction for $\\beta$")

    ax.set_xlabel("$\\beta$")
    ax.set_ylabel("success rate")
    ax.yaxis.set_major_locator(MultipleLocator(0.1))
    ax.grid(axis='y', alpha=0.5)

    save_pdf(fig, path, desc)
    plt.close(fig)


if __name__ == '__main__':
    pass
