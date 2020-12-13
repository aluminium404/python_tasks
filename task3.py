import pylab
import numpy as np
from numpy.fft import fft, fftshift
import matplotlib.pyplot as plt

class GaussianPlaneWave:
    ''' Класс с уравнением плоской волны для гауссова сигнала в дискретном виде
    d - определяет задержку сигнала.
    w - определяет ширину сигнала.
    Sc - число Куранта.
    eps - относительная диэлектрическая проницаемость среды, в которой расположен источник.
    mu - относительная магнитная проницаемость среды, в которой расположен источник.
    '''
    def __init__(self, d, w, Sc=1.0, eps=1.0, mu=1.0):
        self.d = d
        self.w = w
        self.Sc = Sc
        self.eps = eps
        self.mu = mu

    def getE(self, m, q):
        '''
        Расчет поля E в дискретной точке пространства m
        в дискретный момент времени q
        '''
        return np.exp(-(((q - m * np.sqrt(self.eps * self.mu) / self.Sc) - self.d) / self.w) ** 2)

if __name__ == '__main__':
    # Волновое сопротивление свободного пространства
    W0 = 120.0 * np.pi

    # Число Куранта
    Sc = 1.0

    # Скорость света
    c = 3e8

    # Размер ячейки разбиения вдоль оси Х, м
    dx = 1e-3

    # шаг по времени
    dt = Sc * dx / c

    # Параметры гауссова сигнала
    A_0 = 100
    A_max = 100
    F_max = 3e9

    # Ширина и задержка гауссова сигнала
    w_g = np.sqrt(np.log(A_max)) / (np.pi * F_max)
    d_g = w_g * np.sqrt(np.log(A_0))

    
    # Время расчета в отсчетах
    maxTime = 800

    # Размер области моделирования в отсчетах
    maxSize = 500

    # Положение датчика, регистрирующего поле
    probePos = 200

    # Положение источника в отсчетах
    sourcePos = 50

    Ez = np.zeros(maxSize)
    Hy = np.zeros(maxSize)

    source = GaussianPlaneWave(d_g/dt, w_g/dt, Sc)

    # Поле, зарегистрированное в датчике в зависимости от времени
    probeTimeEz = np.zeros(maxTime)
    probeTimeEz[0] = Ez[probePos]

    # Подготовка к отображению поля в пространстве
    xlist = np.arange(0, maxSize*dx, dx)
    
    # Включить интерактивный режим для анимации
    pylab.ion()

    # Создание окна для графика
    fig, ax = pylab.subplots()

    # Установка отображаемых интервалов по осям
    ax.set_xlim(0, maxSize*dx)
    ax.set_ylim(-1.1, 1.1)

    # Установка меток по осям
    ax.set_xlabel('x, м')
    ax.set_ylabel('Ez, В/м')

    # Включить сетку на графике
    ax.grid()

    # Отобразить поле в начальный момент времени
    line, = ax.plot(xlist, Ez)

    # Отобразить положение источника
    ax.plot(sourcePos*dx, 0, 'ok')

    # Отобразить положение датчика
    ax.plot(probePos*dx, 0, 'xr')

    for q in range(1, maxTime):
        # Граничные условия для поля H
        Hy[-1] = Hy[-2]

        # Расчет компоненты поля H
        Ez_shift = Ez[1:]
        Hy[:-1] = Hy[:-1] + (Ez_shift - Ez[:-1]) * Sc / W0

        # Источник возбуждения с использованием метода
        # Total Field / Scattered Field
        Hy[sourcePos - 1] -= (Sc / W0) * source.getE(sourcePos, q)
        # Hy[sourcePos - 1] -= (Sc / W0) * source.getE(0, q)
        
        # Граничные условия для поля E
        Ez[0] = Ez[1]

        # Расчет компоненты поля E
        Hy_shift = Hy[:-1]
        Ez[1:] = Ez[1:] + (Hy[1:] - Hy_shift) * Sc * W0

        # Источник возбуждения с использованием метода
        # Total Field / Scattered Field
        Ez[sourcePos] += Sc * source.getE(sourcePos - 0.5, q + 0.5)
        # Ez[sourcePos] += Sc * source.getE(-0.5, q + 0.5)

        # Регистрация поля в точке
        probeTimeEz[q] = Ez[probePos]

        if q % 2 == 0:
            # Обновить данные на графике
            line.set_ydata(Ez)
            fig.canvas.draw()
            fig.canvas.flush_events()

    # Отключить интерактивный режим по завершению анимации
    pylab.ioff()

    # Расчет спектра
    spectrum = np.abs(fft(probeTimeEz))
    spectrum = fftshift(spectrum)

    # Шаг по частоте
    df = 1.0 / (maxTime * dt)

    # Расчет частоты
    freq = np.arange(-maxTime / 2 * df, maxTime / 2 * df, df)

    #Расчет области
    rang = np.arange(0, maxTime*dt, dt)

    # Отображение сигнала, сохраненного в датчике
    fig, (ax1, ax2) = plt.subplots(2, 1)
    ax1.set_xlim(0, maxTime*dt)
    ax1.set_ylim(-1.1, 1.1)
    ax1.set_xlabel('t, с')
    ax1.set_ylabel('Ez, В/м')
    ax1.plot(rang, probeTimeEz)
    ax1.grid()

    # Отображение спектра
    ax2.plot(freq, spectrum / np.max(spectrum))
    ax2.grid()
    ax2.set_xlabel('Частота, Гц')
    ax2.set_ylabel('|S| / |Smax|')
    ax2.set_xlim(0, maxTime/8*df)
    pylab.show()

