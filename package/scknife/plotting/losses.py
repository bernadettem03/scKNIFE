import matplotlib.pyplot as plt


def plot_losses(loss_history):
    plt.plot(loss_history, marker='o')
    plt.title("Graph-Regularized Sparse NMF Loss")
    plt.xlabel("Iteration")
    plt.ylabel("Objective Value")
    plt.grid(True)
    plt.tight_layout()
    plt.show()
