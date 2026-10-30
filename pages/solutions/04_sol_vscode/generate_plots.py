import numpy as np
import matplotlib.pyplot as plt

# ==========================================
# 1. PARAMETRY WEJŚCIOWE (Dostosuj do indeksu)
# ==========================================
# Przykład dla indeksu 123456:
# Pierwsze dwie cyfry: 12 -> v0 = 12 m/s
# Dwie ostatnie cyfry: 56 -> alpha = 56 deg
v0 = 36.0          # [m/s] Prędkość początkowa (2 pierwsze cyfry indeksu)
alpha_deg = 41.0   # [deg] Kąt wyrzutu (2 ostatnie cyfry indeksu)

g = 9.81           # [m/s^2] Przyspieszenie ziemskie
m = 1.0            # [kg] Masa ciała
b = 0.15           # [kg/s] Współczynnik oporu powietrza (model liniowy)

alpha_rad = np.radians(alpha_deg)

# ==========================================
# 2. RZUT UKOŚNY BEZ OPORU (Analitycznie)
# ==========================================
t_flight_no_drag = (2 * v0 * np.sin(alpha_rad)) / g
t_no_drag = np.linspace(0, t_flight_no_drag, 200)

x_no_drag = v0 * np.cos(alpha_rad) * t_no_drag
y_no_drag = v0 * np.sin(alpha_rad) * t_no_drag - 0.5 * g * t_no_drag**2

vx_no_drag = np.full_like(t_no_drag, v0 * np.cos(alpha_rad))
vy_no_drag = v0 * np.sin(alpha_rad) - g * t_no_drag

# ==========================================
# 3. RZUT UKOŚNY Z OPOREM POWIETRZA (Numerycznie - metoda Eulera)
# ==========================================
dt = 0.001
t_drag = [0.0]
x_drag = [0.0]
y_drag = [0.0]
vx_drag = [v0 * np.cos(alpha_rad)]
vy_drag = [v0 * np.sin(alpha_rad)]

while y_drag[-1] >= 0:
    vx_curr = vx_drag[-1]
    vy_curr = vy_drag[-1]

    # Równania ruchu (opór liniowy F_d = -b*v)
    ax = -(b / m) * vx_curr
    ay = -g - (b / m) * vy_curr

    # Aktualizacja prędkości i położenia
    vx_next = vx_curr + ax * dt
    vy_next = vy_curr + ay * dt
    x_next = x_drag[-1] + vx_curr * dt
    y_next = y_drag[-1] + vy_curr * dt

    x_drag.append(x_next)
    y_drag.append(y_next)
    vx_drag.append(vx_next)
    vy_drag.append(vy_next)
    t_drag.append(t_drag[-1] + dt)

# ==========================================
# 4. GENEROWANIE I ZAPIS WYKRESÓW
# ==========================================
plt.style.use(
    'seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')

# Wykres 1: Trajektoria ruchu y(x)
plt.figure(figsize=(8, 5))
plt.plot(x_no_drag, y_no_drag, 'b-', label='Bez oporu powietrza')
plt.plot(x_drag, y_drag, 'r--', label='Z oporem powietrza')
plt.title('Trajektoria rzutu ukośnego $y(x)$')
plt.xlabel('Położenie $x$ [m]')
plt.ylabel('Wysokość $y$ [m]')
plt.ylim(bottom=0)
plt.legend()
plt.tight_layout()
plt.savefig('trajektoria.png', dpi=300)
plt.close()

# Wykres 2: Wysokość w funkcji czasu y(t)
plt.figure(figsize=(8, 4))
plt.plot(t_no_drag, y_no_drag, 'b-', label='Bez oporu')
plt.plot(t_drag, y_drag, 'r--', label='Z oporem')
plt.title('Wysokość w funkcji czasu $y(t)$')
plt.xlabel('Czas $t$ [s]')
plt.ylabel('Wysokość $y$ [m]')
plt.ylim(bottom=0)
plt.legend()
plt.tight_layout()
plt.savefig('wysokosc_czas.png', dpi=300)
plt.close()

# Wykres 3: Składowe prędkości vx(t) i vy(t)
fig, axs = plt.subplots(1, 2, figsize=(12, 4))

axs[0].plot(t_no_drag, vx_no_drag, 'b-', label='Bez oporu')
axs[0].plot(t_drag, vx_drag, 'r--', label='Z oporem')
axs[0].set_title('Składowa poziom $v_x(t)$')
axs[0].set_xlabel('Czas $t$ [s]')
axs[0].set_ylabel('$v_x$ [m/s]')
axs[0].legend()

axs[1].plot(t_no_drag, vy_no_drag, 'b-', label='Bez oporu')
axs[1].plot(t_drag, vy_drag, 'r--', label='Z oporem')
axs[1].set_title('Składowa pionowa $v_y(t)$')
axs[1].set_xlabel('Czas $t$ [s]')
axs[1].set_ylabel('$v_y$ [m/s]')
axs[1].legend()

plt.tight_layout()
plt.savefig('predkosci.png', dpi=300)
plt.close()

print("Wykresy zostały pomyślnie wygenerowane i zapisane!")
