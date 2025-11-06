"""
1D Axial Model for Packed-Bed LHTES Unit
Replicates Figure 16 from "Modelling a packed-bed latent heat thermal energy storage unit 
and studying its performance using different paraffins" by Klitou et al. (2024)

This model simulates the charging process in a packed-bed LHTES tank using a 1D axial porous media approach 
with local thermal non-equilibrium and lumped PCM assumption. Inlet temperature is coupled via constant heat addition 
from the solar collector to match the integrated system.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

class PackedBedLHTESModel:
    def __init__(self):
        """Initialize the packed-bed LHTES model with parameters from the paper"""
        
        # Geometry parameters
        self.D_tank = 0.36  # Tank diameter [m]
        self.H_tank = 0.47  # Tank height [m]
        self.A = np.pi * (self.D_tank / 2)**2  # Cross-sectional area [m²]
        self.eps = 0.49  # Porosity
        self.d_p = 0.055  # Capsule diameter [m]
        self.nz = 50  # Number of axial nodes
        self.dz = self.H_tank / self.nz  # Axial spacing [m]
        
        # Operating parameters
        self.m_dot = 0.1  # Mass flow rate [kg/s] (6 kg/min for validation)
        self.Q = 375  # Constant heat input from solar collector [W]
        
        # Paraffin properties (matched to validation case)
        self.rho_s = 916  # Solid density [kg/m³]
        self.rho_l = 775  # Liquid density [kg/m³]
        self.cp_s = 2100  # Solid specific heat [J/(kg·K)]
        self.cp_l = 2100  # Liquid specific heat [J/(kg·K)]
        self.k_s = 0.4  # Solid thermal conductivity [W/(m·K)]
        self.k_l = 0.2  # Liquid thermal conductivity [W/(m·K)]
        self.L = 170000  # Latent heat of fusion [J/kg]
        self.Tm = 333.15  # Melting temperature [K] (~60°C)
        self.dT_phase = 2.0  # Phase change temperature range [K]
        
        # Water (HTF) properties
        self.rho_w = 998  # Density [kg/m³]
        self.cp_w = 4182  # Specific heat [J/(kg·K)]
        self.k_w = 0.6  # Thermal conductivity [W/(m·K)]
        self.mu_w = 1e-3  # Dynamic viscosity [Pa·s]
        self.Pr_w = 7  # Prandtl number
        
        # Initial conditions
        self.T_initial = 298.15  # Initial temperature [K] (25°C)
        
        # Calculate flow and heat transfer coefficients
        self.u = self.m_dot / (self.rho_w * self.A)  # Superficial velocity [m/s]
        self.Re = self.rho_w * self.u * self.d_p / self.mu_w  # Reynolds number
        self.Nu = 2 + 1.1 * self.Re**0.6 * self.Pr_w**(1/3)  # Wakao correlation
        self.h_i = self.Nu * self.k_w / self.d_p  # Interstitial HTC [W/(m²·K)]
        self.h_vol = 6 * (1 - self.eps) / self.d_p * self.h_i  # Volumetric HTC [W/(m³·K)]
        
        # Initialize temperature arrays
        self.T_htf = np.ones(self.nz) * self.T_initial
        self.T_pcm = np.ones(self.nz) * self.T_initial
    
    def get_thermal_properties(self, T):
        """Get thermal properties based on temperature (no internal convection enhancement)"""
        if T < (self.Tm - self.dT_phase/2):
            # Solid phase
            rho = self.rho_s
            cp = self.cp_s
            k = self.k_s
            phase = 'solid'
        elif T > (self.Tm + self.dT_phase/2):
            # Liquid phase
            rho = self.rho_l
            cp = self.cp_l
            k = self.k_l
            phase = 'liquid'
        else:
            # Mushy zone
            liquid_fraction = (T - (self.Tm - self.dT_phase/2)) / self.dT_phase
            rho = self.rho_s * (1 - liquid_fraction) + self.rho_l * liquid_fraction
            cp = (self.cp_s + self.cp_l) / 2 + self.L / self.dT_phase
            k = self.k_s * (1 - liquid_fraction) + self.k_l * liquid_fraction
            phase = 'mushy'
        
        return rho, cp, k, phase
    
    def heat_transfer_ode(self, t, y):
        """ODE system for heat transfer in packed-bed LHTES"""
        T_htf = y[:self.nz]
        T_pcm = y[self.nz:]
        
        dTdt_htf = np.zeros(self.nz)
        dTdt_pcm = np.zeros(self.nz)
        
        # Calculate inlet temperature from outlet + collector heat addition
        delta_T = self.Q / (self.m_dot * self.cp_w)
        T_in = T_htf[-1] + delta_T
        
        # HTF equations (upwind scheme for advection)
        for i in range(self.nz):
            if i == 0:
                T_upstream = T_in
            else:
                T_upstream = T_htf[i-1]
            
            # Advection + convection from PCM
            dTdt_htf[i] = - (self.u / self.eps) * (T_htf[i] - T_upstream) / self.dz + \
                          (self.h_vol / (self.rho_w * self.cp_w * self.eps)) * (T_pcm[i] - T_htf[i])
        
        # PCM equations (lumped)
        for i in range(self.nz):
            rho, cp, k, phase = self.get_thermal_properties(T_pcm[i])
            dTdt_pcm[i] = (self.h_vol / (rho * cp * (1 - self.eps))) * (T_htf[i] - T_pcm[i])
        
        return np.append(dTdt_htf, dTdt_pcm)
    
    def simulate(self, t_final=3600*12, dt_save=60):
        """Run the simulation"""
        t_eval = np.arange(0, t_final + dt_save, dt_save)
        
        y0 = np.append(self.T_htf, self.T_pcm)
        
        print("Starting packed-bed LHTES simulation...")
        print(f"Initial temp: {self.T_initial - 273.15:.1f}°C")
        print(f"Melting temp: {self.Tm - 273.15:.1f}°C")
        print(f"Mass flow rate: {self.m_dot * 60:.1f} kg/min")
        print(f"Collector heat input: {self.Q} W")
        print(f"Re_p: {self.Re:.1f}")
        print(f"h_vol: {self.h_vol:.0f} W/(m³·K)")
        
        sol = solve_ivp(self.heat_transfer_ode, [0, t_final], y0, 
                        t_eval=t_eval, method='BDF', rtol=1e-4, atol=1e-6, max_step=10)
        
        if sol.success:
            print("Simulation completed successfully!")
        else:
            print(f"Simulation issues: {sol.message}")
        
        T_htf_history = sol.y[:self.nz, :]
        T_pcm_history = sol.y[self.nz:, :]
        
        return sol.t, T_htf_history, T_pcm_history
    
    def plot_results(self, t, T_htf_history, T_pcm_history):
        """Plot results similar to Figure 16 in the paper"""
        t_hours = t / 3600
        
        # Select three axial positions (bottom, middle, top)
        pos_bottom = 0
        pos_middle = self.nz // 2
        pos_top = self.nz - 1
        
        plt.figure(figsize=(10, 8))
        
        # PCM temperatures
        plt.plot(t_hours, T_pcm_history[pos_bottom, :] - 273.15, 'b-', 
                 label='PCM bottom', linewidth=2)
        plt.plot(t_hours, T_pcm_history[pos_middle, :] - 273.15, 'b--', 
                 label='PCM middle', linewidth=1.5)
        plt.plot(t_hours, T_pcm_history[pos_top, :] - 273.15, 'b:', 
                 label='PCM top', linewidth=1.5)
        
        # HTF temperatures
        plt.plot(t_hours, T_htf_history[pos_bottom, :] - 273.15, 'r-', 
                 label='HTF bottom', linewidth=2)
        plt.plot(t_hours, T_htf_history[pos_middle, :] - 273.15, 'r--', 
                 label='HTF middle', linewidth=1.5)
        plt.plot(t_hours, T_htf_history[pos_top, :] - 273.15, 'r:', 
                 label='HTF top', linewidth=1.5)
        
        # Add melting temperature line
        plt.axhline(y=self.Tm - 273.15, color='k', linestyle='--', alpha=0.3, 
                    label=f'Melting temp ({self.Tm - 273.15:.1f}°C)')
        
        plt.xlabel('Time [h]', fontsize=12)
        plt.ylabel('Temperature [°C]', fontsize=12)
        plt.title('Packed-bed LHTES Temperature Distribution (COMSOL Validation Case)', fontsize=12)
        plt.legend(loc='lower right', fontsize=10)
        plt.grid(True, alpha=0.3)
        plt.xlim(3.5, 6.5)
        plt.ylim(55, 65)
        plt.minorticks_on()
        plt.grid(True, which='minor', alpha=0.1)
        plt.tight_layout()
        plt.savefig('packed_bed_temperature_figure16.png', dpi=150, bbox_inches='tight')
        plt.show()

def main():
    """Main function to run the simulation and generate plot"""
    model = PackedBedLHTESModel()
    
    # Run simulation for 12 hours
    t, T_htf_history, T_pcm_history = model.simulate(t_final=3600*12)
    
    # Plot results
    model.plot_results(t, T_htf_history, T_pcm_history)
    
    # Print summary
    print("\nSimulation Results Summary")
    print(f"Final PCM bottom: {T_pcm_history[0, -1] - 273.15:.1f}°C")
    print(f"Final PCM middle: {T_pcm_history[model.nz//2, -1] - 273.15:.1f}°C")
    print(f"Final PCM top: {T_pcm_history[-1, -1] - 273.15:.1f}°C")
    print(f"Final HTF bottom: {T_htf_history[0, -1] - 273.15:.1f}°C")
    print(f"Final HTF middle: {T_htf_history[model.nz//2, -1] - 273.15:.1f}°C")
    print(f"Final HTF top: {T_htf_history[-1, -1] - 273.15:.1f}°C")

if __name__ == "__main__":
    main()