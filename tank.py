#!python3

"""
TankEvap
---------
An engineering calculation that estimates the evaporative water losses from
storage tanks at a facility. This model uses a modified method established
within Chapter 7 of EPA AP-42.

The method estimates evaporation from both standing and working losses. 

Evaporative Losses From Fixed Roof Tank:
Lt = Ls + Lw
   where:
      Lt = total losses (lb/yr)
      Ls = standing storage losses (lb/yr)
      Lw = working losses (lb/yr)

"""

import math
from pint import UnitRegistry
from thermo.chemical import Chemical

ureg = UnitRegistry()


IDEAL_GAS_CONSTANT = 10.731


class TankEvap:
    """
    The original Ch. 7 calculation.

    Parameters
    ----------
    (reference name) pythonic_name : dtype (eng. units)
    ..........
    (D)   tank_diameter : Float (ft)
    (Hs)  tank_shell_height : Float (ft)
    (Hl)  average_liquid_height : Float (ft)
    (Hx)  maximum_liquid_height : Float (ft)
    (Hro) roof_outage : Float (ft)
    (Mv)  vapor_molecular_weight : Float (lb/lb mol)
    (Tla) daily_avg_liquid_surface_temp : Float (F)
    (a)   tank_paint_solar_absorptance : Float (--)
    (I)   daily_solar_insolation_factor : Float (btu/ft^2-d)
    (Tax) daily_max_ambient_temperature : Float (F)
    (Tan) daily_min_ambient_temperature : Float (F)
    (Pb)  breather_vent_pressure_setting : Float (psia)
    (Q)   net_throughput : Float (gpm)

    Examples
    --------
    Calculating Evaporation Rates

    >>> cst = tank.TankEvap(
            tank_diameter=48,
            tank_shell_height=40,
            average_liquid_height=25.3,
            maximum_liquid_height=36.9,
            roof_outage=0,
            vapor_molecular_weight=18.0152,
            daily_avg_liquid_surface_temp=108,
            tank_paint_solar_absorptance=0.1,
            daily_solar_insolation_factor=1208,
            daily_max_ambient_temperature=63.5,
            daily_min_ambient_temperature=44.5,
            breather_vent_pressure_setting=0,
            net_throughput=520,
            )
    >>> cst.calculate_total_losses()
    >>> cst.total_losses
    32178.064468484252
    """

    def __init__(
        self,
        tank_diameter=None,
        tank_shell_height=None,
        average_liquid_height=None,
        maximum_liquid_height=None,
        roof_outage=None,
        vapor_molecular_weight=None,
        daily_avg_liquid_surface_temp=None,
        tank_paint_solar_absorptance=None,
        daily_solar_insolation_factor=None,
        daily_max_ambient_temperature=None,
        daily_min_ambient_temperature=None,
        breather_vent_pressure_setting=None,
        net_throughput=None,
        atmospheric_pressure=14.7,
    ):
        def degF_to_degR(value):
            Q_ = ureg.Quantity
            return Q_(value, ureg.degF).to("degR").magnitude

        self.tank_diameter = tank_diameter
        self.tank_shell_height = tank_shell_height
        self.average_liquid_height = average_liquid_height
        self.maximum_liquid_height = maximum_liquid_height
        self.roof_outage = roof_outage
        self.vapor_molecular_weight = vapor_molecular_weight
        self.daily_avg_liquid_surface_temp = daily_avg_liquid_surface_temp
        self.tank_paint_solar_absorptance = tank_paint_solar_absorptance
        self.daily_solar_insolation_factor = daily_solar_insolation_factor
        self.daily_max_ambient_temperature = daily_max_ambient_temperature
        self.daily_min_ambient_temperature = daily_min_ambient_temperature
        self.breather_vent_pressure_setting = breather_vent_pressure_setting
        self.net_throughput = net_throughput
        self.daily_avg_liquid_surface_temp_R = degF_to_degR(
            daily_avg_liquid_surface_temp
        )
        self.atmospheric_pressure = atmospheric_pressure

    def calculate_vapor_space_outage(self):
        "Hvo = Hs - Hi + Hro"
        self.vapor_space_outage = (
            self.tank_shell_height - self.average_liquid_height + self.roof_outage
        )

    def calculate_vapor_space_volume(self):
        "Vv = (ⲡ/4)*D^2*Hvo"
        self.calculate_vapor_space_outage()
        self.vapor_space_volume = (
            (math.pi / 4) * self.tank_diameter ** 2 * self.vapor_space_outage
        )

    def calculate_vented_vapor_saturation_factor(self):
        "Ks = 1/(1 + 0.053 * Pva * Hvo)"
        self.calculate_vapor_space_outage()
        self.calculate_vapor_pressures()
        self.vented_vapor_saturation_factor = 1 / (
            1 + 0.053 * self.vapor_pressure_at_tla * self.vapor_space_outage
        )

    def calculate_vapor_density(self):
        "Wv = (Mv * Pva)/(R * Tla)"
        self.calculate_vapor_pressures()
        self.vapor_density = (
            self.vapor_molecular_weight * self.vapor_pressure_at_tla
        ) / (IDEAL_GAS_CONSTANT * self.daily_avg_liquid_surface_temp_R)

    def calculate_daily_vapor_temperature_range(self):
        "dTv = 0.72 * (Tax-Tan) + 0.028 * a * I"
        self.daily_vapor_temperature_range = (
            0.72
            * (self.daily_max_ambient_temperature - self.daily_min_ambient_temperature)
            + 0.028
            * self.tank_paint_solar_absorptance
            * self.daily_solar_insolation_factor
        )

    def calculate_daily_liquid_surface_temp_extremes(self):
        """
        Tlx = Tla + 0.25 * dTv (All Temps in degR)
        Tln = Tla - 0.25 * dTv (All Temps in degR)
        """
        self.calculate_daily_vapor_temperature_range()

        def degR_to_degF(value):
            Q_ = ureg.Quantity
            return Q_(value, ureg.degR).to("degF").magnitude

        tlx = (
            self.daily_avg_liquid_surface_temp_R
            + 0.25 * self.daily_vapor_temperature_range
        )
        self.daily_max_liquid_surface_temp = degR_to_degF(tlx)
        tln = (
            self.daily_avg_liquid_surface_temp_R
            - 0.25 * self.daily_vapor_temperature_range
        )
        self.daily_min_liquid_surface_temp = degR_to_degF(tln)

    def calculate_vapor_pressures(self):
        """
        at Tla: Pva = exp [A - (B/Tla + C)]
        at Tlx: Pvx = exp [A - (B/Tlx + C)]
        at Tln: Pvn = exp [A - (B/Tln + C)]

        Dependency library, thermo, uses degK for vapor press calcs.
        """
        self.calculate_daily_liquid_surface_temp_extremes()

        def degF_to_degK(value):
            Q_ = ureg.Quantity
            return Q_(value, ureg.degF).to("degK").magnitude

        tla = degF_to_degK(self.daily_avg_liquid_surface_temp)
        tlx = degF_to_degK(self.daily_max_liquid_surface_temp)
        tln = degF_to_degK(self.daily_min_liquid_surface_temp)
        self.calculate_daily_vapor_temperature_range()

        def water_vapor_press_at(temp):
            wat = Chemical("water")
            wvp = wat.VaporPressure(temp) * ureg.Pa
            return wvp.to("psi").magnitude

        # Pva
        try:
            self.vapor_pressure_at_tla
        except AttributeError:
            self.vapor_pressure_at_tla = water_vapor_press_at(tla)
        # Pvx
        try:
            self.vapor_pressure_at_tlx
        except AttributeError:
            self.vapor_pressure_at_tlx = water_vapor_press_at(tlx)
        # Pvn
        try:
            self.vapor_pressure_at_tln
        except AttributeError:
            self.vapor_pressure_at_tln = water_vapor_press_at(tln)

    def calculate_vapor_space_expansion_factor(self):
        "Ke = dTv / Tla + (Pvx - Pvn -Pb)/(Pa - Pva)"
        self.calculate_vapor_pressures()
        self.vapor_space_expansion_factor = (
            self.daily_vapor_temperature_range / self.daily_avg_liquid_surface_temp_R
        ) + (
            (
                self.vapor_pressure_at_tlx
                - self.vapor_pressure_at_tln
                - self.breather_vent_pressure_setting
            )
            / (self.atmospheric_pressure - self.vapor_pressure_at_tla)
        )

    def calculate_standing_losses(self):
        "Ls = 365 * Vv * Wv * Ke * Ks"
        self.calculate_vapor_space_volume()
        self.calculate_vented_vapor_saturation_factor()
        self.calculate_vapor_density()
        self.calculate_vapor_space_expansion_factor()
        self.standing_losses = (
            365
            * self.vapor_space_volume
            * self.vapor_density
            * self.vapor_space_expansion_factor
            * self.vented_vapor_saturation_factor
        )

    def calculate_working_losses(self):
        "Lw = 0.001 * Mv * Pva * Q * Kn * Kp"
        self.annual_net_throughput = self.net_throughput * 24 * 60 * 365 / 42
        self.working_tank_volume = (
            math.pi * self.tank_diameter ** 2 * self.maximum_liquid_height / 4
        )
        self.turnovers = 5.615 * self.annual_net_throughput / self.working_tank_volume
        if self.turnovers > 36:
            self.turnover_factor = (180 + self.turnovers) / (6 * self.turnovers)
        else:
            self.turnover_factor = 1
        self.working_loss_product_factor = 1
        self.working_losses = (
            0.001
            * self.vapor_molecular_weight
            * self.vapor_pressure_at_tla
            * self.annual_net_throughput
            * self.turnover_factor
            * self.working_loss_product_factor
        )

    def calculate_total_losses(self):
        "Lt = Ls + Lw"
        self.calculate_standing_losses()
        self.calculate_working_losses()
        self.total_losses = self.standing_losses + self.working_losses


class TankEvap2020(TankEvap):
    """
    The modified Ch. 7 equation from June 2020.

    Parameters
    ----------
    (reference name) pythonic_name : dtype (eng. units)
    ..........
    (D)   tank_diameter : Float (ft)
    (Hs)  tank_shell_height : Float (ft)
    (Hl)  average_liquid_height : Float (ft)
    (Hx)  maximum_liquid_height : Float (ft)
    (Hro) roof_outage : Float (ft)
    (Mv)  vapor_molecular_weight : Float (lb/lb mol)
    (Tla) daily_avg_liquid_surface_temp : Float (F)
    (a)   tank_paint_solar_absorptance : Float (--)
    (I)   daily_solar_insolation_factor : Float (btu/ft^2-d)
    (Tax) daily_max_ambient_temperature : Float (F)
    (Tan) daily_min_ambient_temperature : Float (F)
    (Pb)  breather_vent_pressure_setting : Float (psia)
    (Q)   net_throughput : Float (gpm)

    Examples
    --------
    Calculating Evaporation Rates

    >>> cst = tank.TankEvap2020(
            tank_diameter=48,
            tank_shell_height=40,
            average_liquid_height=25.3,
            maximum_liquid_height=36.9,
            minimum_liquid_height=13.7,
            roof_outage=0,
            vapor_molecular_weight=18.0152,
            daily_avg_liquid_surface_temp=108,
            tank_paint_solar_absorptance=0.1,
            daily_solar_insolation_factor=1208,
            daily_max_ambient_temperature=63.5,
            daily_min_ambient_temperature=44.5,
            breather_vent_pressure_setting=0,
            annual_sum_liquid_increases=0,
            )
    >>> cst.calculate_total_losses()
    >>> cst.total_losses
    32178.064468484252
    """

    def __init__(
        self,
        tank_diameter=None,
        tank_shell_height=None,
        average_liquid_height=None,
        maximum_liquid_height=None,
        minimum_liquid_height=None,
        roof_outage=None,
        vapor_molecular_weight=None,
        daily_avg_liquid_surface_temp=None,
        tank_paint_solar_absorptance=None,
        daily_solar_insolation_factor=None,
        daily_max_ambient_temperature=None,
        daily_min_ambient_temperature=None,
        breather_vent_pressure_setting=None,
        annual_sum_liquid_increases=None,
        atmospheric_pressure=14.7,
    ):
        def degF_to_degR(value):
            Q_ = ureg.Quantity
            return Q_(value, ureg.degF).to("degR").magnitude

        self.tank_diameter = tank_diameter
        self.tank_shell_height = tank_shell_height
        self.average_liquid_height = average_liquid_height
        self.maximum_liquid_height = maximum_liquid_height
        self.minimum_liquid_height = minimum_liquid_height
        self.roof_outage = roof_outage
        self.vapor_molecular_weight = vapor_molecular_weight
        self.daily_avg_liquid_surface_temp = daily_avg_liquid_surface_temp
        self.tank_paint_solar_absorptance = tank_paint_solar_absorptance
        self.daily_solar_insolation_factor = daily_solar_insolation_factor
        self.daily_max_ambient_temperature = daily_max_ambient_temperature
        self.daily_min_ambient_temperature = daily_min_ambient_temperature
        self.breather_vent_pressure_setting = breather_vent_pressure_setting
        self.annual_sum_liquid_increases = annual_sum_liquid_increases
        self.daily_avg_liquid_surface_temp_R = degF_to_degR(
            daily_avg_liquid_surface_temp
        )
        self.atmospheric_pressure = atmospheric_pressure

    def calculate_net_working_loss_throughput(self):
        "Vq = (𝜮Hqi) * (π/4) * D^2"
        self.net_working_loss_throughput = (
            self.annual_sum_liquid_increases * (math.pi / 4) * self.tank_diameter ** 2
        )

    def calculate_working_losses(self):
        "Lw = Vq * Kn * Kp * Wv * Kb"
        self.calculate_vapor_density()
        self.calculate_net_working_loss_throughput()
        self.working_tank_volume = (
            math.pi * self.tank_diameter ** 2 * self.maximum_liquid_height / 4
        )
        self.turnovers = self.net_working_loss_throughput / (
            self.maximum_liquid_height - self.minimum_liquid_height
        )
        if self.turnovers > 36:
            self.turnover_factor = (180 + self.turnovers) / (6 * self.turnovers)
        else:
            self.turnover_factor = 1
        working_loss_product_factor = 1
        vent_setting_correction_factor = 1
        factors = working_loss_product_factor * vent_setting_correction_factor
        self.working_losses = (
            self.net_working_loss_throughput
            * self.turnover_factor
            * self.vapor_density
            * factors
        )
