# Ciało na sprężynie — od prawa Hooke’a do rezonansu

*Raport obliczeniowy i przykład wykorzystania Markdown w dokumentacji naukowej.*

**Temat:** jednowymiarowy ruch masy połączonej ze sprężyną.  
**Metoda:** rozwiązania analityczne, symulacja w Pythonie i interpretacja wykresów.  
**Charakter danych:** przykład modelowy; wszystkie wykresy przedstawiają obliczenia, a nie pomiary laboratoryjne.

Hooke stwierdził, że 
> siła sprężystości jest proporcjonalna do wydłużenia sprężyny. 

Newton stwierdził, że
> przyspieszenie jest proporcjonalne do siły. 

Połączenie tych dwóch praw daje równanie ruchu, którego rozwiązanie opisuje drgania harmoniczne. W tym raporcie pokazujemy, jak wychylenie, prędkość i energia zmieniają się w czasie oraz jak wpływa na nie tłumienie i wymuszenie.

Wychylenie mówi, gdzie znajduje się ciało. Prędkość mówi, co wydarzy się za chwilę. Dopiero razem opisują stan oscylatora.


## Spis treści

1. [Cel i model fizyczny](#1-cel-i-model-fizyczny)
2. [Równanie ruchu i jego rozwiązanie](#2-równanie-ruchu-i-jego-rozwiązanie)
3. [Przykład liczbowy i przebiegi czasowe](#3-przykład-liczbowy-i-przebiegi-czasowe)
4. [Energia i przestrzeń fazowa](#4-energia-i-przestrzeń-fazowa)
5. [Drgania tłumione](#5-drgania-tłumione)
6. [Wymuszenie i rezonans](#6-wymuszenie-i-rezonans)
7. [Kod i odtwarzanie wykresów](#7-kod-i-odtwarzanie-wykresów)
8. [Jak wykonać doświadczenie](#8-jak-wykonać-doświadczenie)
9. [Wnioski](#9-wnioski)
10. [Dalsza lektura](#10-dalsza-lektura)

---

## 1. Cel i model fizyczny

Dlaczego masa po odciągnięciu od położenia równowagi wraca, mija je i ponownie zawraca? Sprężyna wywiera siłę przywracającą, a bezwładność sprawia, że ciało nie zatrzymuje się w równowadze. Powstają **drgania**.

Celem raportu jest opisanie położenia, prędkości, przyspieszenia i energii oraz sprawdzenie, jak zmieniają się one po dodaniu oporu i okresowego wymuszenia.

### 1.1. Założenia

- Ruch odbywa się wzdłuż jednej osi.
- Masa $m$ jest stała, a masę sprężyny pomijamy.
- Sprężyna pracuje w zakresie liniowej sprężystości.
- W modelu idealnym pomijamy opory ruchu.
- W modelu tłumionym przyjmujemy siłę oporu proporcjonalną do prędkości.

Współrzędna $x$ oznacza **wychylenie względem położenia równowagi**, nie całkowitą długość sprężyny. Dla poziomej sprężyny bez dodatkowej stałej siły położenie równowagi odpowiada jej naturalnej długości.

![Schemat poziomego układu: ściana, sprężyna, masa, położenie równowagi i siła skierowana przeciwnie do dodatniego wychylenia](files/sprezyna/schemat.png)

*Rysunek 1. Ciało pokazano dla dodatniego wychylenia. Siła sprężystości działa w lewo. Schemat nie zachowuje skali.*

### 1.2. Symbole i jednostki

| Symbol | Znaczenie | Jednostka SI |
| :--- | :--- | :---: |
| $m$ | masa ciała | $\mathrm{kg}$ |
| $k$ | stała sprężystości | $\mathrm{N/m}$ |
| $x$, $A$ | wychylenie i amplituda | $\mathrm{m}$ |
| $v$, $a$ | prędkość i przyspieszenie | $\mathrm{m/s}$, $\mathrm{m/s^2}$ |
| $\omega_0$ | częstość kołowa drgań własnych | $\mathrm{rad/s}$ |
| $T$, $f$ | okres i częstotliwość | $\mathrm{s}$, $\mathrm{Hz}$ |
| $b$ | współczynnik oporu lepkiego | $\mathrm{N\,s/m}$ |
| $\gamma$, $\zeta$ | parametr tłumienia i względny współczynnik tłumienia | $\mathrm{s^{-1}}$, bez jednostki |
| $F_0$, $\Omega$ | amplituda siły i częstość kołowa wymuszenia | $\mathrm{N}$, $\mathrm{rad/s}$ |

Nie należy mylić $f$ z $\omega_0$: zachodzi $\omega_0=2\pi f$.

## 2. Równanie ruchu i jego rozwiązanie

### 2.1. Prawo Hooke’a i druga zasada Newtona

Siła sprężystości jest proporcjonalna do wychylenia i ma przeciwny zwrot:

$$
F_s=-kx.
$$

Z drugiej zasady Newtona otrzymujemy:

$$
m\frac{d^2x}{dt^2}=-kx.
$$

Po podzieleniu przez masę:

$$
\ddot{x}+\omega_0^2x=0,
\qquad
\omega_0=\sqrt{\frac{k}{m}}.
$$

Kropki oznaczają pochodne po czasie: $\dot{x}=v$ oraz $\ddot{x}=a$. Jest to równanie *oscylatora harmonicznego*. Opis ruchu masy na sprężynie omawia także [OpenStax: Simple Harmonic Motion](https://openstax.org/books/university-physics-volume-1/pages/15-1-simple-harmonic-motion).

### 2.2. Rozwiązanie z warunkami początkowymi

Dla $x(0)=x_0$ oraz $v(0)=v_0$:

$$
x(t)=x_0\cos(\omega_0t)+\frac{v_0}{\omega_0}\sin(\omega_0t).
$$

Równoważnie można zapisać:

$$
x(t)=A\cos(\omega_0t+\varphi),
\qquad
A=\sqrt{x_0^2+\left(\frac{v_0}{\omega_0}\right)^2}.
$$

Dla $A>0$ fazę określa $\varphi=\mathrm{atan2}(-v_0/\omega_0,x_0)$. Funkcja `atan2` uwzględnia właściwą ćwiartkę kąta. Jeśli $A=0$, ciało pozostaje w równowadze i faza jest nieokreślona.

Różniczkowanie daje:

$$
v(t)=-A\omega_0\sin(\omega_0t+\varphi),
$$

$$
a(t)=-A\omega_0^2\cos(\omega_0t+\varphi)=-\omega_0^2x(t).
$$

Stąd maksymalne wartości bezwzględne prędkości i przyspieszenia:

$$
v_{\max}=A\omega_0,
\qquad
a_{\max}=A\omega_0^2.
$$

### 2.3. Okres i wpływ parametrów

$$
T=\frac{2\pi}{\omega_0}=2\pi\sqrt{\frac{m}{k}},
\qquad
f=\frac{1}{T}=\frac{1}{2\pi}\sqrt{\frac{k}{m}}.
$$

W liniowym modelu idealnym okres nie zależy od amplitudy. Czterokrotnie większa masa wydłuża okres dwukrotnie, a czterokrotnie większa stała sprężystości skraca go dwukrotnie.

**Kontrola jednostek:** $k/m$ ma wymiar $\mathrm{s^{-2}}$, więc $\sqrt{k/m}$ ma wymiar $\mathrm{s^{-1}}$.

### 2.4. A jeśli sprężyna jest pionowa?

Niech $y$ będzie wydłużeniem od długości naturalnej, mierzonym w dół. Wtedy:

$$
m\ddot{y}=mg-ky,
\qquad
y_{\mathrm{eq}}=\frac{mg}{k}.
$$

Po wprowadzeniu $x=y-y_{\mathrm{eq}}$ otrzymujemy ponownie:

$$
m\ddot{x}+kx=0.
$$

Stała siła grawitacji przesuwa równowagę. Przy przyjętych założeniach nie zmienia okresu drgań wokół tej równowagi.

---

## 3. Przykład liczbowy i przebiegi czasowe

Przyjmujemy $m=0{,}50\,\mathrm{kg}$, $k=20{,}0\,\mathrm{N/m}$, $x_0=0{,}10\,\mathrm{m}$ oraz $v_0=0$. Ciało zostaje odciągnięte o dziesięć centymetrów i puszczone bez nadania prędkości.

| Wielkość | Wartość | Interpretacja |
| :--- | ---: | :--- |
| Amplituda $A$ | $0{,}100\,\mathrm{m}$ | największe wychylenie |
| Częstość kołowa $\omega_0$ | $6{,}325\,\mathrm{rad/s}$ | tempo zmiany fazy |
| Okres $T$ | $0{,}993\,\mathrm{s}$ | czas pełnego drgania |
| Częstotliwość $f$ | $1{,}007\,\mathrm{Hz}$ | liczba drgań na sekundę |
| Maksymalna szybkość | $0{,}632\,\mathrm{m/s}$ | w położeniu równowagi |
| Maksymalne $\lvert a\rvert$ | $4{,}00\,\mathrm{m/s^2}$ | w skrajnych położeniach |
| Energia całkowita $E$ | $0{,}100\,\mathrm{J}$ | stała w modelu idealnym |

Przebiegi, z $t$ wyrażonym w sekundach, mają postać:

$$
x(t)=0{,}10\cos(\sqrt{40}\,t)\;\mathrm{m},
$$

$$
v(t)=-0{,}10\sqrt{40}\sin(\sqrt{40}\,t)\;\mathrm{m/s},
$$

$$
a(t)=-4{,}00\cos(\sqrt{40}\,t)\;\mathrm{m/s^2}.
$$

![Trzy wykresy w funkcji czasu: wychylenie, prędkość i przyspieszenie ciała na sprężynie w ciągu trzech okresów](files/sprezyna/polozenie_predkosc_przyspieszenie.png)

*Rysunek 2. Wspólna oś czasu ułatwia porównanie faz. Szare linie oznaczają kolejne pełne okresy.*

Podczas jednego okresu:

1. W chwili $t=0$ ciało znajduje się w $x=A$, ma $v=0$ i przyspiesza w lewo.
2. W chwili $t=T/4$ przechodzi przez $x=0$ z prędkością $v=-A\omega_0$.
3. W chwili $t=T/2$ dociera do $x=-A$, zatrzymuje się chwilowo i przyspiesza w prawo.
4. W chwili $t=3T/4$ wraca przez równowagę z $v=+A\omega_0$.
5. W chwili $t=T$ powraca do początkowego stanu.

~~W położeniu równowagi ciało musi się zatrzymać.~~ W tym przykładzie właśnie tam porusza się najszybciej; zerowe jest przyspieszenie.

## 4. Energia i przestrzeń fazowa

### 4.1. Zamiana energii

Energia potencjalna sprężystości oraz energia kinetyczna wynoszą:

$$
E_p=\frac{1}{2}kx^2,
\qquad
E_k=\frac{1}{2}mv^2.
$$

Dla $x=A\cos(\omega_0t+\varphi)$:

$$
E_p(t)=\frac{1}{2}kA^2\cos^2(\omega_0t+\varphi),
$$

$$
E_k(t)=\frac{1}{2}kA^2\sin^2(\omega_0t+\varphi).
$$

Ponieważ $\sin^2\theta+\cos^2\theta=1$:

$$
E=E_p+E_k=\frac{1}{2}kA^2=\mathrm{const}.
$$

Można sprawdzić to bez znajomości rozwiązania:

$$
\frac{dE}{dt}=mv\dot{v}+kx\dot{x}=v(m\ddot{x}+kx)=0.
$$

Średnie po pełnym okresie spełniają $\langle E_k\rangle=\langle E_p\rangle=E/2$. Każdy rodzaj energii zmienia się z okresem $T/2$, ponieważ zależy od kwadratu sinusa lub cosinusa.

![Energia potencjalna i kinetyczna wymieniają się okresowo, a ich suma pozostaje równa 0,1 J](files/sprezyna/energia.png)

*Rysunek 3. Energia całkowita jest stała. W skrajnych położeniach cała energia jest potencjalna; w równowadze — kinetyczna.*

### 4.2. Portret fazowy

**Przestrzeń fazowa** opisuje stan układu parą $(x,v)$. W mechanice Hamiltona zwykle używa się pary $(x,p)$, gdzie $p=mv$; dla stałej masy jest to jedynie przeskalowanie osi pionowej.

Eliminując czas z rozwiązania, otrzymujemy:

$$
\frac{x^2}{A^2}+\frac{v^2}{A^2\omega_0^2}=1.
$$

Jest to elipsa, której półosie wynoszą $A$ i $A\omega_0$. Po normalizacji $X=x/A$ oraz $V=v/(A\omega_0)$ równanie przyjmuje postać $X^2+V^2=1$, czyli opisuje okrąg.

![Portret fazowy: elipsa w osiach wychylenie–prędkość i okrąg po normalizacji, ze strzałkami wskazującymi kierunek ruchu](files/sprezyna/przestrzen_fazowa.png)

*Rysunek 4. Jeden zamknięty obieg odpowiada jednemu okresowi. Punkt początkowy znajduje się po prawej stronie; ruch przebiega zgodnie z ruchem wskazówek zegara.*

Z portretu można odczytać:

- **zamkniętą orbitę** — układ wraca do tego samego stanu;
- **przecięcia z osią $x$** — punkty zwrotne, w których $v=0$;
- **przecięcia z osią $v$** — przejścia przez równowagę;
- **rozmiar orbity** — większa amplituda oznacza większą energię.

### 4.3. Zapis macierzowy

Równanie drugiego rzędu można zastąpić układem dwóch równań pierwszego rzędu:

$$
\dot{x}=v,
\qquad
\dot{v}=-\omega_0^2x.
$$

W zapisie macierzowym:

$$
\frac{d}{dt}
\begin{pmatrix}
x \\
v \\
\end{pmatrix}
= \begin{pmatrix}
0 & 1 \\
-\omega_0^2 & 0 \\
\end{pmatrix}
\begin{pmatrix}
x \\
v \\
\end{pmatrix}.
$$

Wartości własne macierzy to $\lambda=\pm i\omega_0$. Ich czysto urojona postać odpowiada oscylacjom bez zaniku. Taki zapis jest też naturalnym punktem wyjścia do całkowania numerycznego.

---

## 5. Drgania tłumione

### 5.1. Opór proporcjonalny do prędkości

Przyjmujemy $F_{\mathrm{op}}=-bv$, gdzie $b>0$. Równanie przyjmuje postać:

$$
m\ddot{x}+b\dot{x}+kx=0.
$$

Wprowadzamy parametry:

$$
\gamma=\frac{b}{2m},
\qquad
\zeta=\frac{b}{2\sqrt{mk}}=\frac{\gamma}{\omega_0}.
$$

Wtedy $\ddot{x}+2\gamma\dot{x}+\omega_0^2x=0$. Rozwiązania klasyfikuje się według $\zeta$, co opisuje [OpenStax: Damped Oscillations](https://openstax.org/books/university-physics-volume-1/pages/15-5-damped-oscillations).

### 5.2. Słabe tłumienie: $0<\zeta<1$

Częstość kołowa drgań tłumionych jest mniejsza od własnej:

$$
\omega_d=\sqrt{\omega_0^2-\gamma^2},
\qquad
T_d=\frac{2\pi}{\omega_d}.
$$

Dla dowolnych warunków początkowych:

$$
x(t)=e^{-\gamma t}
\left[
x_0\cos(\omega_dt)
+\frac{v_0+\gamma x_0}{\omega_d}\sin(\omega_dt)
\right].
$$

Można też zapisać $x(t)=A_d e^{-\gamma t}\cos(\omega_dt+\varphi_d)$, gdzie:

$$
A_d=\sqrt{x_0^2+\left(\frac{v_0+\gamma x_0}{\omega_d}\right)^2}.
$$

Obwiednie wynoszą $\pm A_d e^{-\gamma t}$. Nawet gdy $v_0=0$, zazwyczaj $A_d\ne x_0$: składnik sinusowy jest potrzebny, aby początkowa prędkość rzeczywiście była zerowa.

Dla $b=0{,}40\,\mathrm{N\,s/m}$ otrzymujemy $\gamma=0{,}40\,\mathrm{s^{-1}}$, $\zeta\approx0{,}0632$, $\omega_d\approx6{,}312\,\mathrm{rad/s}$ oraz $T_d\approx0{,}995\,\mathrm{s}$.

![Zanikające wychylenie z wykładniczymi obwiedniami oraz spirala w przestrzeni fazowej zmierzająca do równowagi](files/sprezyna/tlumienie.png)

*Rysunek 5. Układ traci energię: wychylenie maleje, a orbita fazowa zwija się do punktu $(0,0)$.*

### 5.3. Gdzie znika energia?

$$
\frac{dE}{dt}=v(m\ddot{x}+kx)=-bv^2\leq0.
$$

Energia mechaniczna przechodzi w energię wewnętrzną otoczenia. Jej chwilowa wartość nie jest ogólnie dokładnie równa $E(0)e^{-2\gamma t}$. Dla słabego tłumienia uśredniony przebieg można przybliżyć zanikiem wykładniczym:

$$
\overline{E}(t)\approx\overline{E}(0)e^{-2\gamma t}.
$$

Czas spadku obwiedni amplitudy do $1/e$ wartości początkowej wynosi $\tau_A=1/\gamma$, a analogiczny czas dla uśrednionej energii to około $\tau_E=1/(2\gamma)$. W tym przykładzie są to odpowiednio $2{,}50\,\mathrm{s}$ i $1{,}25\,\mathrm{s}$.

### 5.4. Trzy rodzaje tłumienia

| Zakres | Nazwa | Postać rozwiązania ogólnego | Zachowanie |
| :--- | :--- | :--- | :--- |
| $0<\zeta<1$ | podkrytyczne | zanikający sinus i cosinus | oscylacje wokół równowagi |
| $\zeta=1$ | krytyczne | $(C_1+C_2t)e^{-\omega_0t}$ | brak oscylacji |
| $\zeta>1$ | nadkrytyczne | $C_1e^{r_1t}+C_2e^{r_2t}$ | brak oscylacji; wolna składowa zaniku |

Dla tłumienia nadkrytycznego:

$$
r_{1,2}=-\gamma\pm\sqrt{\gamma^2-\omega_0^2}.
$$

Stałe $C_1$ i $C_2$ wynikają z warunków początkowych. Dla puszczenia z $x_0>0$ i $v_0=0$ rozwiązanie krytyczne ma szczególnie prostą postać:

$$
x(t)=x_0(1+\omega_0t)e^{-\omega_0t}.
$$

![Porównanie wychylenia dla tłumienia podkrytycznego, krytycznego i nadkrytycznego przy tych samych warunkach początkowych](files/sprezyna/rodzaje_tlumienia.png)

*Rysunek 6. Dla pokazanych warunków początkowych tłumienie krytyczne zapewnia szybszy powrót bez oscylacji niż nadkrytyczne. Zwiększanie oporu nie zawsze przyspiesza powrót.*

## 6. Wymuszenie i rezonans

### 6.1. Okresowa siła zewnętrzna

Jeżeli układ napędza siła $F(t)=F_0\cos(\Omega t)$:

$$
m\ddot{x}+b\dot{x}+kx=F_0\cos(\Omega t).
$$

Dla $b>0$ rozwiązanie składa się z zanikającej części przejściowej i trwałej odpowiedzi okresowej. Po zaniku części przejściowej:

$$
x_{\mathrm{ust}}(t)=A_{\mathrm{ust}}(\Omega)\cos(\Omega t-\delta),
$$

$$
A_{\mathrm{ust}}(\Omega)=
\frac{F_0}{\sqrt{(k-m\Omega^2)^2+(b\Omega)^2}},
$$

$$
\delta=\mathrm{atan2}(b\Omega,k-m\Omega^2).
$$

Opóźnienie fazowe $\delta$ rośnie od $0$ do $\pi$, a dla $\Omega=\omega_0$ wynosi $\pi/2$. Przy małej częstości amplituda zbliża się do statycznego wychylenia $F_0/k$; przy dużej maleje jak $F_0/(m\Omega^2)$.

Podstawowe zjawisko rezonansu omawia [OpenStax: Forced Oscillations](https://openstax.org/books/university-physics-volume-1/pages/15-6-forced-oscillations). Poniżej rozróżniamy maksimum **amplitudy wychylenia** od innych miar odpowiedzi.

### 6.2. Częstość rezonansowa i dobroć

Maksimum amplitudy wychylenia występuje przy dodatniej częstości:

$$
\Omega_r=\sqrt{\omega_0^2-2\gamma^2}
=\omega_0\sqrt{1-2\zeta^2},
\qquad
0<\zeta<\frac{1}{\sqrt{2}}.
$$

Gdy $\zeta\geq1/\sqrt{2}$, amplituda nie ma maksimum dla dodatniej częstości i maleje od wartości statycznej. Przy słabym tłumieniu $\Omega_r\approx\omega_0$, ale te wielkości nie są dokładnie równe.

Dobroć oscylatora w tym modelu definiujemy jako:

$$
Q=\frac{m\omega_0}{b}=\frac{1}{2\zeta}.
$$

Dla $\Omega=\omega_0$:

$$
A_{\mathrm{ust}}(\omega_0)=\frac{F_0}{b\omega_0}
=Q\frac{F_0}{k}.
$$

![Charakterystyki rezonansowe: amplituda wychylenia i opóźnienie fazowe w funkcji względnej częstości wymuszenia dla trzech wartości tłumienia](files/sprezyna/rezonans.png)

*Rysunek 7. Przyjęto $F_0=0{,}20\,\mathrm{N}$. Słabsze tłumienie daje wyższy i węższy pik amplitudy. Pionowa linia oznacza $\Omega/\omega_0=1$, a nie dokładne położenie każdego maksimum.*

### 6.3. Moc i granica bez tłumienia

W stanie ustalonym średnia moc dostarczana przez siłę zewnętrzną równoważy straty:

$$
\langle P\rangle=\frac{1}{2}b\Omega^2A_{\mathrm{ust}}^2.
$$

Dla stałej amplitudy siły i $b>0$ maksimum tej mocy występuje przy $\Omega=\omega_0$. Jest to inna charakterystyka niż amplituda wychylenia, której maksimum przypada na $\Omega_r$.

Przy $b=0$ i dokładnym rezonansie nie istnieje ograniczona odpowiedź ustalona. Dla $x(0)=v(0)=0$:

$$
x(t)=\frac{F_0}{2m\omega_0}t\sin(\omega_0t).
$$

Obwiednia rośnie liniowo w czasie. W rzeczywistym układzie taki wzrost ograniczają między innymi opory i nieliniowość sprężyny.

---

## 7. Kod i odtwarzanie wykresów

Wykresy wykonano za pomocą NumPy, Matplotlib i SciPy. Wszystkie osie mają opisy i jednostki, a kolory są spójne w całym raporcie. Pliki PNG zapisujemy w rozdzielczości 180 dpi.

### 7.1. Uruchomienie w Colabie lub lokalnie

1. Otwórz [Google Colab](https://colab.research.google.com/) i utwórz notebook.
2. Skopiuj **cztery poniższe bloki `python`** do czterech kolejnych komórek kodowych.
3. Uruchom je po kolei od bloku A do D. Kolejne bloki korzystają ze zmiennych i funkcji z wcześniejszych.
4. Pobierz katalog `files/sprezyna` z wynikami i umieść go obok tego raportu, zachowując strukturę podfolderów.

Lokalnie możesz uruchomić dołączony [skrypt `generuj_wykresy.py`](generuj_wykresy.py), zawierający te same cztery bloki. Z katalogu `pages/examples/01_markdown` wykonaj:

```bash
python -m pip install numpy matplotlib scipy
python generuj_wykresy.py
```

Program zapisuje obrazy względem bieżącego katalogu roboczego. W Colabie będzie to zwykle `/content`; lokalnie uruchamiaj go z katalogu raportu. Funkcja `zapisz` wyświetla każdy wykres także w notebooku.

### 7.2. Blok A — parametry, styl i schemat

```python
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.integrate import solve_ivp

OUT = Path("files/sprezyna")
OUT.mkdir(parents=True, exist_ok=True)
BLUE, ORANGE, GREEN = "#2463A6", "#D47722", "#168578"
PURPLE, INK = "#8455A4", "#233044"
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 11,
    "axes.titlesize": 13, "axes.titleweight": "bold",
    "axes.labelcolor": INK, "text.color": INK,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.alpha": 0.18,
    "lines.linewidth": 2.3, "figure.facecolor": "white",
    "savefig.facecolor": "white",
})

m, k = 0.50, 20.0
x0, v0 = 0.10, 0.0
w0 = np.sqrt(k / m)
T = 2 * np.pi / w0
A = np.hypot(x0, v0 / w0)

def zapisz(fig, nazwa):
    fig.savefig(OUT / f"{nazwa}.png", dpi=180)
    if plt.get_backend().lower() != "agg":
        plt.show()
    plt.close(fig)

print(f"omega_0 = {w0:.6f} rad/s; T = {T:.6f} s")
print(f"f = {1/T:.6f} Hz; E = {0.5*k*A**2:.6f} J")

fig, ax = plt.subplots(figsize=(10, 3.4), layout="constrained")
ax.set(xlim=(-0.5, 7.0), ylim=(-0.6, 2.2))
ax.axis("off")
ax.add_patch(Rectangle((0, 0.3), 0.22, 1.25,
                       facecolor="#CBD5E1", hatch="///", edgecolor=INK))
ax.plot([0.2, 6.6], [0.3, 0.3], color=INK, lw=1.4)
sx = np.linspace(0.55, 4.25, 500)
ax.plot([0.22, 0.55], [0.95, 0.95], color=BLUE)
ax.plot(sx, 0.95 + 0.20*np.sin(np.linspace(0, 18*np.pi, 500)), color=BLUE)
ax.plot([4.25, 4.6], [0.95, 0.95], color=BLUE)
ax.add_patch(Rectangle((4.6, 0.3), 1.1, 1.3,
                       facecolor="#E3EEF9", edgecolor=BLUE, lw=2))
ax.text(5.15, 0.95, "$m$", ha="center", va="center", fontsize=22)
ax.text(2.3, 1.45, "sprężyna o stałej $k$", ha="center")
ax.plot([3.9, 3.9], [0.1, 1.65], "--", color="#64748B", lw=1.2)
ax.annotate("", xy=(5.15, -0.13), xytext=(3.9, -0.13),
            arrowprops={"arrowstyle": "<->", "color": INK})
ax.text(4.52, -0.42, "$x>0$", ha="center")
ax.text(3.9, 1.85, "równowaga: $x=0$", ha="center")
ax.annotate("", xy=(4.0, 1.7), xytext=(5.7, 1.7),
            arrowprops={"arrowstyle": "->", "color": ORANGE, "lw": 2.5})
ax.text(5.3, 1.98, "$F_s=-kx$", color=ORANGE, ha="center")
ax.set_title("Model: masa na poziomej sprężynie", loc="left")
zapisz(fig, "schemat")
```

### 7.3. Blok B — ruch idealny, energia i portret fazowy

```python
t = np.linspace(0, 3*T, 1801)
x = x0*np.cos(w0*t) + (v0/w0)*np.sin(w0*t)
v = -x0*w0*np.sin(w0*t) + v0*np.cos(w0*t)
a = -w0**2*x
Ep, Ek = 0.5*k*x**2, 0.5*m*v**2

fig, axes = plt.subplots(3, 1, figsize=(10, 7), sharex=True,
                         layout="constrained")
for ax, y, label, color in zip(
    axes, [x, v, a], ["x [m]", "v [m/s]", "a [m/s²]"],
    [BLUE, ORANGE, PURPLE]
):
    ax.plot(t, y, color=color)
    ax.set_ylabel(label)
    ax.axhline(0, color=INK, lw=0.7, alpha=0.5)
    for n in range(1, 4):
        ax.axvline(n*T, color=INK, ls=":", lw=0.9, alpha=0.4)
axes[0].set_title("Ruch idealny: trzy okresy drgań", loc="left")
axes[-1].set_xlabel("Czas t [s]")
zapisz(fig, "polozenie_predkosc_przyspieszenie")

fig, ax = plt.subplots(figsize=(10, 4.2), layout="constrained")
ax.plot(t, Ep, color=BLUE, label="Potencjalna")
ax.plot(t, Ek, color=ORANGE, label="Kinetyczna")
ax.plot(t, Ep+Ek, color=GREEN, ls="--", label="Całkowita")
ax.set(xlabel="Czas t [s]", ylabel="Energia [J]", ylim=(-0.005, 0.12))
ax.set_title("Energia zmienia postać, ale jej suma pozostaje stała", loc="left")
ax.legend(loc="upper right", ncol=3, fontsize=9)
zapisz(fig, "energia")

# Jeden okres wystarczy, aby narysować pełną orbitę.
tf = np.linspace(0, T, 601)
xf = x0*np.cos(w0*tf) + (v0/w0)*np.sin(w0*tf)
vf = -x0*w0*np.sin(w0*tf) + v0*np.cos(w0*tf)
fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), layout="constrained")
for ax, xx, yy in zip(axes, [xf, xf/A], [vf, vf/(A*w0)]):
    ax.plot(xx, yy, color=BLUE)
    ax.scatter(xx[0], yy[0], color=ORANGE, zorder=3, label="Start")
    for j in [40, 190, 340, 490]:
        ax.annotate("", xy=(xx[j+15], yy[j+15]), xytext=(xx[j], yy[j]),
                    arrowprops={"arrowstyle": "->", "color": BLUE, "lw": 2})
    ax.axhline(0, color=INK, lw=0.7, alpha=0.4)
    ax.axvline(0, color=INK, lw=0.7, alpha=0.4)
axes[0].set(xlabel="Wychylenie x [m]", ylabel="Prędkość v [m/s]",
            title="Portret fazowy w jednostkach SI", xticks=np.linspace(-A, A, 5))
axes[0].legend(loc="upper right")
axes[1].set(xlabel="x / A", ylabel="v / (Aω₀)", title="Po normalizacji: okrąg")
axes[1].set_aspect("equal", adjustable="box")
zapisz(fig, "przestrzen_fazowa")

# Sprawdzamy prawo zachowania energii, a nie tylko wygląd rysunku.
blad_E = np.max(np.abs(Ep+Ek - 0.5*k*A**2))
assert blad_E < 1e-12
print(f"Maksymalny bezwzględny błąd energii: {blad_E:.3e} J")
```

### 7.4. Blok C — tłumienie i rozwiązanie numeryczne

Metoda `solve_ivp` oblicza rozwiązanie układu $\dot{x}=v$, $\dot{v}=-(b/m)v-(k/m)x$. Jej argumenty opisano w [dokumentacji SciPy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html). Zmienna `t_eval` określa chwile zapisu wyniku; solver sam dobiera wewnętrzne kroki całkowania.

```python
b = 0.40
gamma = b / (2*m)
wd = np.sqrt(w0**2 - gamma**2)
td = np.linspace(0, 8*T, 3201)
C = x0
D = (v0 + gamma*x0) / wd
q = C*np.cos(wd*td) + D*np.sin(wd*td)
xd = np.exp(-gamma*td)*q
vd = np.exp(-gamma*td)*(
    -gamma*q - C*wd*np.sin(wd*td) + D*wd*np.cos(wd*td)
)
obwiednia = np.hypot(C, D)*np.exp(-gamma*td)

fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), layout="constrained")
axes[0].plot(td, xd, color=BLUE, label="x(t)")
axes[0].plot(td, obwiednia, "--", color=ORANGE, label="Obwiednie")
axes[0].plot(td, -obwiednia, "--", color=ORANGE)
axes[0].set(xlabel="Czas t [s]", ylabel="Wychylenie x [m]",
            title="Zanikające drgania")
axes[0].legend()
axes[1].plot(xd, vd, color=PURPLE)
axes[1].scatter(xd[0], vd[0], color=ORANGE, zorder=3, label="Start")
axes[1].scatter(0, 0, color=INK, marker="+", s=90, label="Równowaga")
for j in [90, 480, 1000]:
    axes[1].annotate("", xy=(xd[j+22], vd[j+22]), xytext=(xd[j], vd[j]),
                     arrowprops={"arrowstyle": "->", "color": PURPLE, "lw": 2})
axes[1].set(xlabel="Wychylenie x [m]", ylabel="Prędkość v [m/s]",
            title="Spirala w przestrzeni fazowej", xticks=np.linspace(-A, A, 5))
axes[1].legend(fontsize=9)
zapisz(fig, "tlumienie")

def symuluj(zeta, czasy):
    b_test = 2*zeta*np.sqrt(m*k)
    def rhs(czas, stan):
        xx, vv = stan
        return [vv, -(b_test/m)*vv - (k/m)*xx]
    sol = solve_ivp(rhs, (czasy[0], czasy[-1]), [x0, v0],
                    t_eval=czasy, rtol=1e-9, atol=1e-11)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.y

tr = np.linspace(0, 2.5*T, 1201)
fig, ax = plt.subplots(figsize=(10, 4.4), layout="constrained")
for zeta, color, nazwa in [
    (0.15, BLUE, "Podkrytyczne"),
    (1.0, GREEN, "Krytyczne"),
    (2.0, ORANGE, "Nadkrytyczne"),
]:
    xr, vr = symuluj(zeta, tr)
    ax.plot(tr, xr, color=color, label=f"{nazwa}: ζ = {zeta:g}")
ax.axhline(0, color=INK, lw=0.8)
ax.set(xlabel="Czas t [s]", ylabel="Wychylenie x [m]")
ax.set_title("Ten sam stan początkowy, trzy rodzaje tłumienia", loc="left")
ax.legend(loc="upper right", fontsize=10)
zapisz(fig, "rodzaje_tlumienia")

# Niezależna kontrola: solver numeryczny kontra wzór analityczny.
xn, vn = symuluj(gamma/w0, td)
blad_x = np.max(np.abs(xn-xd))
blad_v = np.max(np.abs(vn-vd))
assert blad_x < 1e-8 and blad_v < 1e-7
Ed = 0.5*k*xd**2 + 0.5*m*vd**2
assert np.max(np.diff(Ed)) < 1e-12
print(f"Błąd x: {blad_x:.3e} m; błąd v: {blad_v:.3e} m/s")
print("Energia tłumionego oscylatora nie rośnie.")
```

### 7.5. Blok D — charakterystyka rezonansu

```python
F0 = 0.20
r = np.linspace(0, 2.2, 2201)
Omega = r*w0
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), layout="constrained")
for zeta, color in [(0.05, BLUE), (0.15, ORANGE), (0.40, GREEN)]:
    bz = 2*zeta*np.sqrt(m*k)
    amplituda = F0/np.hypot(k-m*Omega**2, bz*Omega)
    faza = np.arctan2(bz*Omega, k-m*Omega**2)
    axes[0].plot(r, 100*amplituda, color=color, label=f"ζ = {zeta:.2f}")
    axes[1].plot(r, faza, color=color)
    # Znane analityczne maksimum powinno zgadzać się z siatką.
    r_max = r[np.argmax(amplituda)]
    assert abs(r_max-np.sqrt(1-2*zeta**2)) < 2*(r[1]-r[0])
for ax in axes:
    ax.axvline(1, color=INK, lw=1, ls="--", alpha=0.5)
    ax.set_xlabel("Względna częstość wymuszenia Ω / ω₀")
    ax.set_xlim(0, 2.2)
axes[0].set(ylabel="Amplituda wychylenia [cm]", title="Odpowiedź amplitudowa")
axes[0].legend()
axes[1].set(ylabel="Opóźnienie fazowe δ [rad]", title="Odpowiedź fazowa",
            yticks=[0, np.pi/2, np.pi], yticklabels=["0", "π/2", "π"])
zapisz(fig, "rezonans")
print(f"Gotowe: zapisano 7 rysunków w {OUT.resolve()}")
```

### 7.6. Struktura plików

```text
pages/examples/01_markdown/
├── raport_cialo_na_sprezynie.md
├── generuj_wykresy.py
└── files/
    └── sprezyna/
        ├── schemat.png
        ├── polozenie_predkosc_przyspieszenie.png
        ├── energia.png
        ├── przestrzen_fazowa.png
        ├── tlumienie.png
        ├── rodzaje_tlumienia.png
        └── rezonans.png
```

**Ważne dla GitHuba:** obrazy są osadzone za pomocą ścieżek względnych, np. `files/sprezyna/energia.png`. Dodając raport do repozytorium, trzeba dodać również pliki graficzne. Notebook Colaba nie odczytuje automatycznie takich ścieżek z repozytorium; obrazy trzeba pobrać albo ponownie wygenerować.

## 8. Jak wykonać doświadczenie

### 8.1. Prosta procedura pomiarowa

1. Zmierz masę ciała $m$.
2. Dla sprężyny pionowej zmierz jej wydłużenie równowagowe $\Delta l$ względem długości naturalnej.
3. Oszacuj stałą sprężystości z zależności $k=mg/\Delta l$.
4. Odciągnij ciało nieznacznie od równowagi i puść bez nadania prędkości.
5. Zmierz czas $t_N$ obejmujący $N$ pełnych okresów, np. $N=20$.
6. Oblicz okres $T_{\mathrm{pom}}=t_N/N$ i porównaj go z $T_{\mathrm{teor}}=2\pi\sqrt{m/k}$.

Względną różnicę możesz zapisać jako:

$$
\varepsilon_T=
\frac{\lvert T_{\mathrm{pom}}-T_{\mathrm{teor}}\rvert}{T_{\mathrm{teor}}}\cdot100\%.
$$

Pomiar wielu okresów ogranicza względny wpływ błędu reakcji przy uruchamianiu i zatrzymywaniu stopera. Warto powtórzyć serię kilka razy, a nie oceniać zgodności na podstawie jednego pomiaru.

### 8.2. Niepewność i wyznaczanie parametrów

Jeśli niepewność pomiaru czasu to $u(t_N)$, a liczbę okresów znamy dokładnie:

$$
u(T_{\mathrm{pom}})=\frac{u(t_N)}{N}.
$$

Dla niezależnych niepewności masy i stałej sprężystości, przy liniowym przybliżeniu propagacji:

$$
\frac{u(T_{\mathrm{teor}})}{T_{\mathrm{teor}}}
\approx\frac{1}{2}
\sqrt{\left(\frac{u(m)}{m}\right)^2+\left(\frac{u(k)}{k}\right)^2}.
$$

Jeśli $k$ wyznaczono z tej samej masy przez $k=mg/\Delta l$, wielkości $m$ i $k$ są skorelowane. Wówczas wygodniej użyć bezpośrednio $T=2\pi\sqrt{\Delta l/g}$ i propagować niepewność tego wzoru.

Z pomiarów dla różnych mas można wyznaczyć $k$, dopasowując prostą:

$$
T^2=\frac{4\pi^2}{k}m.
$$

Z kolei amplitudy kolejnych dodatnich maksimów $a_n$ i $a_{n+1}$ w ruchu podkrytycznym pozwalają oszacować tłumienie:

$$
\Lambda=\ln\left(\frac{a_n}{a_{n+1}}\right)=\gamma T_d,
\qquad
b=2m\frac{\Lambda}{T_d}.
$$

Tutaj $\Lambda$ jest dekrementem logarytmicznym. Korzystamy z maksimów oddzielonych pełnym okresem, nie z sąsiednich ekstremów o przeciwnych znakach.

### 8.3. Ograniczenia modelu

- **Sprężyna rzeczywista:** duże wydłużenie może wykraczać poza zakres prawa Hooke’a.
- **Masa sprężyny:** gdy nie jest mała w porównaniu z $m$, prosty model wymaga korekty.
- **Opory:** tarcie suche nie jest opisane siłą $-bv$; zanikanie amplitudy może mieć inny przebieg.
- **Dodatkowe ruchy:** kołysanie boczne i skręcanie utrudniają opis jednowymiarowy.
- **Rezonans:** bardzo duże amplitudy na wykresie liniowego modelu nie gwarantują jego poprawności dla rzeczywistej sprężyny.

## 9. Wnioski

Układ masa–sprężyna ma częstość własną $\omega_0=\sqrt{k/m}$. W ruchu idealnym wychylenie zmienia się sinusoidalnie, energia całkowita pozostaje stała, a stan układu zatacza zamkniętą orbitę w przestrzeni fazowej.

Tłumienie odbiera energię z szybkością $bv^2$. Przy tłumieniu podkrytycznym amplituda zanika wykładniczo, a portret fazowy ma postać spirali. Tłumienie krytyczne i nadkrytyczne eliminują oscylacje, ale większy opór nie musi oznaczać szybszego powrotu do równowagi.

Siła okresowa podtrzymuje drgania. Dla słabego tłumienia rezonans prowadzi do dużej amplitudy, a jej maksimum leży nieco poniżej częstości własnej. **Maksimum amplitudy wychylenia i maksimum średniej mocy nie występują ogólnie przy tej samej częstości.**

### Lista kontrolna raportu

- [x] Zdefiniowano model, współrzędne i jednostki.
- [x] Wyprowadzono rozwiązanie oraz wzory na energię.
- [x] Dołączono przebiegi czasowe i portrety fazowe.
- [x] Porównano tłumienie i charakterystyki rezonansowe.
- [x] Udostępniono kod pozwalający odtworzyć rysunki.
- [x] Kod sprawdza zachowanie energii i zgodność rozwiązania analitycznego z numerycznym.
- [ ] Raport przeczytany przez czytelnika.

---

## 10. Dalsza lektura

- [Wikipedia: Oscylator harmoniczny](https://pl.wikipedia.org/wiki/Oscylator_harmoniczny) — krótki przegląd pojęcia i zastosowań.
- [OpenStax: Simple Harmonic Motion](https://openstax.org/books/university-physics-volume-1/pages/15-1-simple-harmonic-motion) — podstawowy model drgań.
- [OpenStax: Damped Oscillations](https://openstax.org/books/university-physics-volume-1/pages/15-5-damped-oscillations) — ruch z oporami.
- [OpenStax: Forced Oscillations](https://openstax.org/books/university-physics-volume-1/pages/15-6-forced-oscillations) — wymuszenie i rezonans.

## Nota od autora

Plik jest przykładem połączenia opisu w Markdown, obliczeń w Pythonie oraz wyników zapisanych w repozytorium. Czytelnik może zmienić parametry w bloku A, uruchomić kod ponownie i porównać nowy wynik z przewidywaniami wzorów.
