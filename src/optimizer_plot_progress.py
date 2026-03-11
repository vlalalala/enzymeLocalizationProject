import optuna
import matplotlib.pyplot as plt
import numpy as np

def plot_optimization_progress(storage_path, study_name, product_to_maximize):
    study = optuna.load_study(
        study_name=study_name,
        storage=f"sqlite:///{storage_path}"
    )

    # group trial values by round
    rounds = {}
    for trial in study.trials:
        if trial.value is None:  # skip incomplete trials
            continue
        round_idx = trial.user_attrs.get("round", trial.number // N_TRIALS)  # fallback to inference
        rounds.setdefault(round_idx, []).append(trial.value)

    fig, ax = plt.subplots(figsize=(10, 5))

    best_so_far = []
    current_best = -np.inf

    for round_idx in sorted(rounds.keys()):
        values = rounds[round_idx]
        ax.scatter(
            [round_idx] * len(values),
            values,
            color="steelblue",
            alpha=0.6,
            zorder=2
        )
        current_best = max(current_best, max(values))
        best_so_far.append((round_idx, current_best))

    # best-so-far line across rounds
    best_rounds, best_values = zip(*best_so_far)
    ax.step(best_rounds, best_values, where="post",
            color="red", linewidth=2, label="Best so far")

    ax.set_xlabel("Round")
    ax.set_ylabel(product_to_maximize)
    ax.set_title(f"Optimization progress — {product_to_maximize}")
    ax.set_xticks(sorted(rounds.keys()))
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("optimization_progress.png", dpi=150)
    plt.show()