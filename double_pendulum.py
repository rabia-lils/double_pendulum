import numpy as np
import matplotlib.pyplot as plt

# Fiziksel sabitler
g = 9.81
m1 = 1.0
m2 = 1.0
L1 = 1.0
L2 = 1.0

# Zaman parametreleri
dt = 0.01
T = 20
t = np.arange(0, T, dt)

# Diziler
theta1 = np.zeros(len(t))
theta2 = np.zeros(len(t))
omega1 = np.zeros(len(t))
omega2 = np.zeros(len(t))

# Başlangıç koşulları (radyan)
theta1[0] = np.pi / 2
theta2[0] = np.pi / 2 + 0.01  # küçük fark = kaos
omega1[0] = 0
omega2[0] = 0

# Hareket denklemleri (Euler yöntemi)
for i in range(len(t)-1):
    delta = theta2[i] - theta1[i]

    den1 = (m1 + m2) * L1 - m2 * L1 * np.cos(delta)**2
    den2 = (L2 / L1) * den1

    a1 = (
        m2 * L1 * omega1[i]**2 * np.sin(delta) * np.cos(delta)
        + m2 * g * np.sin(theta2[i]) * np.cos(delta)
        + m2 * L2 * omega2[i]**2 * np.sin(delta)
        - (m1 + m2) * g * np.sin(theta1[i])
    ) / den1

    a2 = (
        -m2 * L2 * omega2[i]**2 * np.sin(delta) * np.cos(delta)
        + (m1 + m2) * g * np.sin(theta1[i]) * np.cos(delta)
        - (m1 + m2) * L1 * omega1[i]**2 * np.sin(delta)
        - (m1 + m2) * g * np.sin(theta2[i])
    ) / den2

    omega1[i+1] = omega1[i] + a1 * dt
    omega2[i+1] = omega2[i] + a2 * dt
    theta1[i+1] = theta1[i] + omega1[i] * dt
    theta2[i+1] = theta2[i] + omega2[i] * dt

# Kartesyen koordinatlar
x1 = L1 * np.sin(theta1)
y1 = -L1 * np.cos(theta1)
x2 = x1 + L2 * np.sin(theta2)
y2 = y1 - L2 * np.cos(theta2)

# Grafik
plt.plot(x2, y2)
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.title("Çift Sarkaç – Kaotik Hareket")
plt.grid()

plt.savefig("double_pendulum_trajectory.png")
plt.show()
