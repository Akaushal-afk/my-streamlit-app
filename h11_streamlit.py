import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

st.title("H11 Graph Calculator")

# Input fields
a = st.number_input("Radius of tread ring (m)", format="%.6f")
b = st.number_input("Tread width (m)", format="%.6f")
rho_A = st.number_input("Density of tread ring (kg/m)", format="%.6f")
m_wheel = st.number_input("Mass of wheel (kg)", format="%.6f")
I_r = st.number_input("Moment of inertia of wheel (kg·m²)", format="%.6f")
EI = st.number_input("Flexural rigidity of tread ring (N·m²)", format="%.6f")
k_r = st.number_input("Radial fundamental spring rate (N/m²)", format="%.6f")
k_t = st.number_input("Circumferential fundamental spring rate (N/m²)", format="%.6f")
p_0 = st.number_input("Inflation pressure (Pa)", format="%.6f")
sigma_theta_0_A = st.number_input("Circumferential tension of tread ring (N)", format="%.6f")
zeta = st.number_input("Damping coefficient", format="%.6f")
Omega = st.number_input("Angular velocity (rad/s)", format="%.6f")

def t_m0_11(s):
    m_0 = rho_A 
    k_0R = k_t * a
    k_R = k_t * a**2
    m_r = I_r / (2 * np.pi * a)
    c_0 = 2 * zeta * np.sqrt(m_0 * k_0R)
    c_r = 2 * zeta * np.sqrt(m_r * k_R)
    num = m_r * s**2 + c_r * s + k_R
    den = (m_0 * s**2 + c_0 * s + k_0R) * (m_r * s**2 + c_r * s + k_R) - k_0R * 2
    return num / den

def t_m1_11(s):
    m_1 = rho_A * 2
    m_a = m_wheel / (np.pi * a)
    k_1 = k_r + k_t - rho_A * Omega**2 * 2
    k_a = k_r + k_t - (m_wheel * Omega**2) / (np.pi * a)
    c_1 = 2 * zeta * np.sqrt(m_1 * k_1)
    c_a = 2 * zeta * np.sqrt(m_a * k_1)
    num = m_a * s**2 + c_a * s + k_a
    den = (m_1 * s**2 + c_1 * s + k_1) * (m_a * s**2 + c_a * s + k_a) - k_1**2
    return num / den

def t_mn_11(n, s):
    if n == 0:
        return t_m0_11(s)
    elif n == 1:
        return t_m1_11(s)
    else:
        m_n = rho_A * (1 + n**2)
        k_n = ((EI / a**4) * n**2 + (sigma_theta_0_A / a**2)) * (n**2 - 1)**2 + (p_0 * b / a) * (n**2 - 1) + k_r * n**2 + k_t - rho_A * Omega**2 * (n**2 + 1)
        c_n = 2 * zeta * np.sqrt(m_n * k_n)
        g_n = 2 * rho_A * n * Omega * (n**2 - 1)
        num = m_n * s**2 + c_n * s + k_n
        den = (m_n * s**2 + c_n * s + k_n)**2 + (g_n * s)**2
        return num / den

def h11(omega, n_max=10):
    s = 1j * omega
    return (0.5 * t_m0_11(s) + sum(t_mn_11(n, s) for n in range(1, n_max + 1))) / (np.pi * a)

if st.button("Graph It!!"):
    try:
        frequencies = np.linspace(0.1, 300, 200)
        h11_magnitudes = [abs(h11(2 * np.pi * f)) for f in frequencies]

        fig, ax = plt.subplots()
        ax.loglog(frequencies, h11_magnitudes)
        ax.set_title('H11 Transfer Function Magnitude')
        ax.set_xlabel('Frequency (Hz)')
        ax.set_ylabel('|H11| Magnitude')
        ax.grid(True, which="both", ls="--")

        st.pyplot(fig)

    except Exception as e:
        st.error(f"Error in calculation: {e}")

if st.button("Clear All"):
    # Streamlit does not have a native clear all button functionality for inputs.
    # Users can manually refresh the page or we can provide instructions.
    st.info("To clear inputs, please refresh the page.")
