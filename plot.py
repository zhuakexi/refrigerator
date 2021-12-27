import matplotlib.pyplot as plt
from matplotlib import colors
def plot_ps_curve(cvd_agg):
    f, ax = plt.subplots(1,1)
    ax.loglog(
        cvd_agg['s_bp'],
        cvd_agg['count.avg'],
    )
    ax.set(
        xlabel='separation, bp',
        ylabel='IC contact frequency')
    ax.set_aspect(1.0)
    ax.grid(lw=0.5)
def plot_lb_ps_curve(lb_cvd_agg):
    f, ax = plt.subplots(1,1)
    ax.loglog(
        lb_cvd_agg['s_bp'],
        lb_cvd_agg['count.avg'],
    )
    ax.set(
        xlabel='separation, bp',
        ylabel='IC contact frequency')
    ax.set_aspect(1.0)
    ax.grid(lw=0.5)