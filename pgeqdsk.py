#!/usr/bin/env python

from dataclasses import dataclass
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as spline
from scipy.interpolate import interp1d
from scipy.interpolate import RectBivariateSpline
import os


@dataclass(eq=False)
class pgeqdsk(object):
    cocos: np.number
    nw: np.number
    nh: np.number
    rdim: np.number
    zdim: np.number
    rcentr: np.number
    rleft: np.number
    zmid: np.number
    rmaxis: np.number
    zmaxis: np.number
    simag: np.number
    sibry: np.number
    bcentr: np.number
    current: np.number
    fpol: np.ndarray
    pres: np.ndarray
    ffprim: np.ndarray
    pprime: np.ndarray
    psirz: np.ndarray
    qpsi: np.ndarray
    nbbbs: np.number
    limitr: np.number
    rbbbs: np.ndarray
    zbbbs: np.ndarray
    rlim: np.ndarray
    zlim: np.ndarray

    def set_cocos(self, cocos):
        if self.cocos not in [3, 11]:
            raise ValueError(f"COCOS conversion is available only from 3 and 11")
        if cocos not in [3, 11]:
            raise ValueError(f"COCOS conversion is available only to 3 and 11")

        if (self.cocos == 3) and (cocos == 11):
            self.simag *= -2.0 * np.pi
            self.sibry *= -2.0 * np.pi
            self.ffprim /= -2.0 * np.pi
            self.pprime /= -2.0 * np.pi
            self.psirz *= -2.0 * np.pi
            self.qpsi = np.abs(self.qpsi) * np.sign(self.current) * np.sign(self.bcentr)

        elif (self.cocos == 11) and (cocos == 3):
            self.simag /= -2.0 * np.pi
            self.sibry /= -2.0 * np.pi
            self.ffprim *= -2.0 * np.pi
            self.pprime *= -2.0 * np.pi
            self.psirz /= -2.0 * np.pi
            self.qpsi = np.abs(self.qpsi)

        self.cocos = cocos

    @property
    def r(self):
        return np.linspace(self.rleft, self.rleft + self.rdim, self.nw)

    @property
    def z(self):
        return np.linspace(
            self.zmid - 0.5 * self.zdim, self.zmid + 0.5 * self.zdim, self.nh
        )

    @property
    def psin(self):
        return np.linspace(0.0, 1.0, self.nw)

    @property
    def psi(self):
        return self.psin * (self.sibry - self.simag) + self.simag

    @property
    def phi(self):
        if self.cocos < 10:
            return np.array(
                [
                    spline(
                        self.psin, 2.0 * np.pi * self.qpsi * (self.sibry - self.simag)
                    ).integral(0.0, p)
                    for p in self.psin
                ]
            )
        else:
            return np.array(
                [
                    spline(self.psin, self.qpsi * (self.sibry - self.simag)).integral(
                        0.0, p
                    )
                    for p in self.psin
                ]
            )

    @property
    def rho(self):
        return np.sqrt(self.phi / (np.pi * self.bcentr))

    @property
    def rhon(self):
        return self.rho / self.rho[-1]

    @property
    def psirzn(self):
        return (self.psirz - self.simag) / (self.sibry - self.simag)

    @property
    def rout(self):
        psinz = interp1d(self.z, self.psirzn)(self.zmid)
        outboard = self.r >= self.rmaxis
        rout = spline(psinz[outboard], self.r[outboard])(self.psin)
        rout[rout < self.rmaxis] = self.rmaxis
        return rout

    @property
    def r_outboard(self):
        psinz = interp1d(self.z, self.psirzn)(self.zmaxis)
        outboard = self.r >= self.rmaxis
        rout = spline(psinz[outboard], self.r[outboard])(self.psin)
        rout[rout < self.rmaxis] = self.rmaxis
        return rout

    def write(self, filename):
        from fgeqdsk import fgeqdsk

        fgeqdsk.nw = self.nw
        fgeqdsk.nh = self.nh
        fgeqdsk.rdim = self.rdim
        fgeqdsk.zdim = self.zdim
        fgeqdsk.rcentr = self.rcentr
        fgeqdsk.rleft = self.rleft
        fgeqdsk.zmid = self.zmid
        fgeqdsk.rmaxis = self.rmaxis
        fgeqdsk.zmaxis = self.zmaxis
        fgeqdsk.simag = self.simag
        fgeqdsk.sibry = self.sibry
        fgeqdsk.bcentr = self.bcentr
        fgeqdsk.current = self.current
        fgeqdsk.fpol = self.fpol
        fgeqdsk.pres = self.pres
        fgeqdsk.ffprim = self.ffprim
        fgeqdsk.pprime = self.pprime
        fgeqdsk.psirz = self.psirz
        fgeqdsk.qpsi = self.qpsi
        fgeqdsk.nbbbs = self.nbbbs
        fgeqdsk.limitr = self.limitr
        fgeqdsk.rbbbs = self.rbbbs
        fgeqdsk.zbbbs = self.zbbbs
        fgeqdsk.rlim = self.rlim
        fgeqdsk.zlim = self.zlim

        fgeqdsk.write_geqdsk(filename)

    @staticmethod
    def xy2rt(x, y):
        r = np.linalg.norm([x, y], axis=0)
        t = np.arctan2(y, x)
        ind = np.argsort(t)
        return r[ind], t[ind]

    @staticmethod
    def from_file(filename, cocos=3):
        if not os.path.isfile(filename):
            raise FileNotFoundError(f"No such file or directory: '{filename}'")

        from fgeqdsk import fgeqdsk

        fgeqdsk.read_geqdsk(filename)
        return pgeqdsk(
            **{
                "cocos": cocos,
                "nw": np.int32(fgeqdsk.nw),
                "nh": np.int32(fgeqdsk.nh),
                "rdim": np.float64(fgeqdsk.rdim),
                "zdim": np.float64(fgeqdsk.zdim),
                "rcentr": np.float64(fgeqdsk.rcentr),
                "rleft": np.float64(fgeqdsk.rleft),
                "zmid": np.float64(fgeqdsk.zmid),
                "rmaxis": np.float64(fgeqdsk.rmaxis),
                "zmaxis": np.float64(fgeqdsk.zmaxis),
                "simag": np.float64(fgeqdsk.simag),
                "sibry": np.float64(fgeqdsk.sibry),
                "bcentr": np.float64(fgeqdsk.bcentr),
                "current": np.float64(fgeqdsk.current),
                "fpol": np.array(fgeqdsk.fpol, dtype=np.float64),
                "pres": np.array(fgeqdsk.pres, dtype=np.float64),
                "ffprim": np.array(fgeqdsk.ffprim, dtype=np.float64),
                "pprime": np.array(fgeqdsk.pprime, dtype=np.float64),
                "psirz": np.array(fgeqdsk.psirz, dtype=np.float64, order="C"),
                "qpsi": np.array(fgeqdsk.qpsi, dtype=np.float64),
                "nbbbs": np.int32(fgeqdsk.nbbbs),
                "limitr": np.int32(fgeqdsk.limitr),
                "rbbbs": np.array(fgeqdsk.rbbbs, dtype=np.float64),
                "zbbbs": np.array(fgeqdsk.zbbbs, dtype=np.float64),
                "rlim": np.array(fgeqdsk.rlim, dtype=np.float64),
                "zlim": np.array(fgeqdsk.zlim, dtype=np.float64),
            }
        )


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import sys

    try:
        filename = sys.argv[1]
    except IndexError:
        print(f"Usage: {sys.argv[0]} $GFILE")
        exit()

    g = pgeqdsk.from_file(filename)
    plt.plot(g.rmaxis, g.zmaxis, "kx")
    plt.contour(g.rdim, g.zdim, g.psirz.T)
    plt.axis("scaled")
    plt.show()
