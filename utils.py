from matplotlib import cm
from matplotlib.collections import LineCollection
import matplotlib.pyplot as plt
import numpy as np
from constants import *
from icecream import ic


def transformProj2Orig(sref, xref, yref, psiref, si, ni, alpha, v):
    tracklength = sref[-1]
    si = si % tracklength
    idxmindist = findClosestS(si, sref)
    idxmindist2 = findSecondClosestS(si, sref, idxmindist)
    t = (si - sref[idxmindist]) / (sref[idxmindist2] - sref[idxmindist])
    x0 = (1 - t) * xref[idxmindist] + t * xref[idxmindist2]
    y0 = (1 - t) * yref[idxmindist] + t * yref[idxmindist2]
    psi0 = (1 - t) * psiref[idxmindist] + t * psiref[idxmindist2]

    x = x0 - ni * np.sin(psi0)
    y = y0 + ni * np.cos(psi0)
    psi = psi0 + alpha
    v = v
    return x, y, psi, v


def findClosestS(si, sref):
    # Get number of elements
    if np.isscalar(si):
        N = 1
    else:
        N = np.array(si).shape[0]
    mindist = 100000 * np.ones(N)
    idxmindist = np.zeros(N)
    for i in range(sref.size):
        di = abs(si - sref[i])
        idxmindist = np.where(di < mindist, i, idxmindist)
        mindist = np.where(di < mindist, di, mindist)
    idxmindist = np.where(idxmindist == sref.size, 1, idxmindist)
    idxmindist = np.where(idxmindist < 1, sref.size - 1, idxmindist)
    return idxmindist.astype(int)


def findSecondClosestS(si, sref, idxmindist):
    d1 = abs(si - sref[idxmindist - 1])  # distance to node before
    d2 = abs(si - sref[(idxmindist + 1) % sref.size])  # distance to node after
    idxmindist2 = np.where(
        d1 > d2, idxmindist + 1, idxmindist - 1
    )  # decide which node is closer
    idxmindist2 = np.where(
        idxmindist2 == sref.size, 0, idxmindist2
    )  # if chosen node is too large
    idxmindist2 = np.where(
        idxmindist2 < 0, sref.size - 1, idxmindist2
    )  # if chosen node is too small

    return idxmindist2


def transformOrig2Proj(sref, xref, yref, psiref, x, y, psi, v):
    idxmindist = findClosestPoint(x, y, xref, yref)
    idxmindist2 = findClosestNeighbour(x, y, xref, yref, idxmindist)
    t = findProjection(x, y, xref, yref, sref, idxmindist, idxmindist2)
    s0 = (1 - t) * sref[idxmindist] + t * sref[idxmindist2]
    x0 = (1 - t) * xref[idxmindist] + t * xref[idxmindist2]
    y0 = (1 - t) * yref[idxmindist] + t * yref[idxmindist2]
    psi0 = (1 - t) * psiref[idxmindist] + t * psiref[idxmindist2]

    s = s0
    n = np.cos(psi0) * (y - y0) - np.sin(psi0) * (x - x0)
    alpha = psi - psi0
    v = v
    return s, n, alpha, v


def findProjection(x, y, xref, yref, sref, idxmindist, idxmindist2):
    vabs = abs(sref[idxmindist] - sref[idxmindist2])
    vl = np.empty(2)
    u = np.empty(2)
    vl[0] = xref[idxmindist2] - xref[idxmindist]
    vl[1] = yref[idxmindist2] - yref[idxmindist]
    u[0] = x - xref[idxmindist]
    u[1] = y - yref[idxmindist]
    t = (vl[0] * u[0] + vl[1] * u[1]) / vabs / vabs
    return t


def findClosestPoint(x, y, xref, yref):
    mindist = 1
    idxmindist = 0
    for i in range(xref.size):
        dist = dist2D(x, xref[i], y, yref[i])
        if dist < mindist:
            mindist = dist
            idxmindist = i
    return idxmindist


def findClosestNeighbour(x, y, xref, yref, idxmindist):
    distBefore = dist2D(x, xref[idxmindist - 1], y, yref[idxmindist - 1])
    distAfter = dist2D(x, xref[idxmindist + 1], y, yref[idxmindist + 1])
    if distBefore < distAfter:
        idxmindist2 = idxmindist - 1
    else:
        idxmindist2 = idxmindist + 1
    if idxmindist2 < 0:
        idxmindist2 = xref.size - 1
    elif idxmindist == xref.size:
        idxmindist2 = 0
    return idxmindist2


def dist2D(x1, x2, y1, y2):
    return np.sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2))


def plot_traj(
    s_ref,
    X_ref,
    Y_ref,
    psi_ref,
    simX,
):
    s = simX[:, 0]
    n = simX[:, 1]
    alpha = simX[:, 2]
    v = simX[:, 3]
    distance = 0.12
    # transform data
    [x, y, _, _] = transformProj2Orig(s_ref, X_ref, Y_ref, psi_ref, s, n, alpha, v)

    # Setup plot
    plt.figure()
    plt.ylim(bottom=-1.75, top=0.35)
    plt.xlim(left=-1.1, right=1.6)
    plt.ylabel(r"$Y$ [m]")
    plt.xlabel(r"$X$ [m]")

    # Plot center line
    plt.plot(X_ref, Y_ref, "--", color="k")

    # Draw Trackboundaries
    Xboundleft = X_ref - distance * np.sin(psi_ref)
    Yboundleft = Y_ref + distance * np.cos(psi_ref)
    Xboundright = X_ref + distance * np.sin(psi_ref)
    Yboundright = Y_ref - distance * np.cos(psi_ref)
    plt.plot(Xboundleft, Yboundleft, color="k", linewidth=1)
    plt.plot(Xboundright, Yboundright, color="k", linewidth=1)

    # Draw driven trajectory
    norm = plt.Normalize(v.min(), v.max())
    points = np.array([x, y]).transpose().reshape(-1, 1, 2)
    segs = np.concatenate([points[:-1], points[1:]], axis=1)
    lc = LineCollection(segs, cmap="brg", norm=norm)
    lc.set_array(v)
    lc.set_linewidth(2)
    line = plt.gca().add_collection(lc)
    plt.gcf().colorbar(line, ax=plt.gca(), label=r"$v$ [m/s]")
    plt.axis("equal")
    plt.title("Driven trajectory")

    # Put markers for s values
    # xi = np.zeros(9)
    # yi = np.zeros(9)
    # xi1 = np.zeros(9)
    # yi1 = np.zeros(9)
    # xi2 = np.zeros(9)
    # yi2 = np.zeros(9)
    # for i in range(int(s_ref[-1]) + 1):
    #     try:
    #         k = list(s_ref).index(i + min(abs(s_ref - i)))
    #     except:
    #         k = list(s_ref).index(i - min(abs(s_ref - i)))
    #     [_, nrefi, _, _] = transformOrig2Proj(
    #         s_ref, X_ref, Y_ref, psi_ref, X_ref[k], Y_ref[k], psi_ref[k], 0
    #     )
    #     [xi[i], yi[i], _, _] = transformProj2Orig(
    #         s_ref, X_ref, Y_ref, psi_ref, s_ref[k], nrefi + 0.24, 0, 0
    #     )
    #     # plt.text(xi[i], yi[i], f'{i}m', fontsize=12,horizontalalignment='center',verticalalignment='center')
    #     plt.text(
    #         xi[i],
    #         yi[i],
    #         "{}m".format(i),
    #         fontsize=12,
    #         horizontalalignment="center",
    #         verticalalignment="center",
    #     )
    #     [xi1[i], yi1[i], _, _] = transformProj2Orig(
    #         s_ref, X_ref, Y_ref, psi_ref, s_ref[k], nrefi + 0.12, 0, 0
    #     )
    #     [xi2[i], yi2[i], _, _] = transformProj2Orig(
    #         s_ref, X_ref, Y_ref, psi_ref, s_ref[k], nrefi + 0.15, 0, 0
    #     )
    #     plt.plot([xi1[i], xi2[i]], [yi1[i], yi2[i]], color="black")
    plt.tight_layout()


def plot_cl(simX, simU, t):
    # plot results
    fig, axs = plt.subplots(2, 1, sharex=True)
    axs[0].step(t, simU[:, 0], color="r", where="post")
    axs[0].step(t, simU[:, 1], color="g", where="post")
    axs[0].legend(
        [r"$dT$", r"$d\delta$"]
        if constant_controls
        else [r"$\dot{T}$", r"$\dot{\delta}$"]
    )
    axs[0].set_ylabel("u")
    axs[0].grid(True)
    if constant_controls:
        axs[1].plot(t, simX[:, :-nu])
        axs[1].step(t, simX[:, -nu:], where="post")
    else:
        axs[1].plot(t, simX[:, :])
    axs[1].set_ylabel("x")
    axs[1].set_xlabel("t")
    axs[1].legend([r"$s$", r"$n$", r"$\alpha$", r"$v$", r"$T$", r"$\delta$"])
    axs[1].grid(True)
    fig.suptitle("Closed-loop simulation")
    fig.tight_layout()


def plot_accels(t, simX, simU, alat_fun, along_fun):
    Nsim = t.shape[0]
    fig, axs = plt.subplots(2, 1, sharex=True)
    alat = np.zeros(Nsim)
    along = np.zeros(Nsim)
    for i in range(Nsim):
        alat[i] = alat_fun(simX[i, :], simU[i, :])
        along[i] = along_fun(simX[i, :], simU[i, :])

    axs[0].plot(t, alat)
    axs[0].plot([t[0], t[-1]], [alat_min, alat_min], "k--")
    axs[0].plot([t[0], t[-1]], [alat_max, alat_max], "k--")
    axs[0].legend([r"$a_{\perp}$", r"$a_\perp$ min/max"])
    axs[0].set_ylabel(r"$a_{\perp}$ [m/s^2]")

    axs[1].plot(t, along)
    axs[1].plot([t[0], t[-1]], [along_min, along_min], "k--")
    axs[1].plot([t[0], t[-1]], [along_max, along_max], "k--")
    axs[1].legend([r"$a_{\parallel}$", r"$a_\parallel$ min/max"])
    axs[1].set_xlabel("t")
    axs[1].set_ylabel(r"$a_{\parallel}$ [m/s^2]")

    fig.suptitle("Accelerations")
    fig.tight_layout()
