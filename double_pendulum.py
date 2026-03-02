import numpy as np
import matplotlib.pyplot as plt

# Fiziksel sabitler
g = 9.81
m1 = 1.0
m2 = 1.0
L1 = 1.0
L2 = 1.0

# Zaman parametreleri
dt = 0.001 
T = 20
t = np.arange(0, T, dt)

# Diziler
theta1 = np.zeros(len(t))
theta2 = np.zeros(len(t))
omega1 = np.zeros(len(t))
omega2 = np.zeros(len(t))

# Başlangıç koşulları
theta1[0] = np.pi / 2
theta2[0] = np.pi / 2 + 0.01 
omega1[0] = 0
omega2[0] = 0

# Hareket denklemleri
for i in range(len(t)-1):
    delta = theta2[i] - theta1[i]

    # Payda hesaplamaları
    den1 = (m1 + m2) * L1 - m2 * L1 * np.cos(delta)**2
    den2 = (L2 / L1) * den1

    # İvmeler (a1 = d(omega1)/dt, a2 = d(omega2)/dt)
    a1 = (m2 * g * np.sin(theta2[i]) * np.cos(delta) 
          - m2 * L1 * omega1[i]**2 * np.sin(delta) * np.cos(delta)
          - m2 * L2 * omega2[i]**2 * np.sin(delta)
          - (m1 + m2) * g * np.sin(theta1[i])) / den1

    a2 = ((m1 + m2) * (L1 * omega1[i]**2 * np.sin(delta) 
          - g * np.sin(theta2[i]) 
          + g * np.sin(theta1[i]) * np.cos(delta)) 
          + m2 * L2 * omega2[i]**2 * np.sin(delta) * np.cos(delta)) / den2
    omega1[i+1] = omega1[i] + a1 * dt
    omega2[i+1] = omega2[i] + a2 * dt
    theta1[i+1] = theta1[i] + omega1[i+1] * dt
    theta2[i+1] = theta2[i] + omega2[i+1] * dt

# Kartesyen koordinat
x1 = L1 * np.sin(theta1)
y1 = -L1 * np.cos(theta1)
x2 = x1 + L2 * np.sin(theta2)
y2 = y1 - L2 * np.cos(theta2)

# Grafik
plt.figure(figsize=(8, 8))
plt.plot(x2, y2, lw=0.5, color='royalblue')
plt.plot(0, 0, 'ko') # Sabit nokta
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.title("Çift Sarkaç - Alt Ucun İzlediği Kaotik Yol")
plt.grid(alpha=0.3)
plt.axis('equal')
plt.show()
