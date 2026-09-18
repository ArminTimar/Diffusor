"""Citation registry.

Every equation, coefficient, constant, buffer and convention in Diffusor
carries a *citation key* that must exist in :data:`REFERENCES`.  A test
(``tests/test_references.py``) enforces this, and ``scripts/build_references.py``
renders ``references.bib`` and ``REFERENCES.md`` from this single source.

Use :func:`cite` inside docstrings/notes to build a compact "Author (year), eq. N"
string, and :func:`format_reference` for a full bibliographic entry.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Optional


@dataclass(frozen=True)
class Reference:
    key: str
    authors: str
    year: int
    title: str
    journal: str = ""
    volume: str = ""
    pages: str = ""
    doi: str = ""
    note: str = ""
    kind: str = "article"   # article | book | incollection | misc

    def short(self) -> str:
        first = self.authors.split(",")[0].split(" and ")[0].strip()
        n_auth = self.authors.count(" and ") + 1
        if n_auth == 1:
            who = first
        elif n_auth == 2:
            second = self.authors.split(" and ")[1].split(",")[0].strip()
            who = f"{first} & {second}"
        else:
            who = f"{first} et al."
        return f"{who} ({self.year})"

    def full(self) -> str:
        s = f"{self.authors} ({self.year}) {self.title}."
        if self.journal:
            s += f" {self.journal}"
            if self.volume:
                s += f" {self.volume}"
            if self.pages:
                s += f":{self.pages}"
            s += "."
        if self.doi:
            s += f" https://doi.org/{self.doi}"
        if self.note:
            s += f" [{self.note}]"
        return s


def _r(**kw) -> Reference:
    return Reference(**kw)


REFERENCES: Dict[str, Reference] = {r.key: r for r in [
    # ---- textbooks / reviews ------------------------------------------------
    _r(key="crank1975", authors="Crank, J.", year=1975,
       title="The Mathematics of Diffusion, 2nd edition", journal="Oxford University Press, Oxford",
       pages="414 pp", kind="book"),
    _r(key="costa2008", authors="Costa, F. and Dohmen, R. and Chakraborty, S.", year=2008,
       title="Time scales of magmatic processes from modelling the zoning patterns of crystals",
       journal="Reviews in Mineralogy and Geochemistry", volume="69", pages="545-594",
       doi="10.2138/rmg.2008.69.14"),
    _r(key="costa2020", authors="Costa, F. and Shea, T. and Ubide, T.", year=2020,
       title="Diffusion chronometry and the timescales of magmatic processes",
       journal="Nature Reviews Earth & Environment", volume="1", pages="201-214",
       doi="10.1038/s43017-020-0038-x"),
    _r(key="dohmen2017", authors="Dohmen, R. and Faak, K. and Blundy, J. D.", year=2017,
       title="Chronometry and speedometry of magmatic processes using chemical diffusion in olivine, plagioclase and pyroxenes",
       journal="Reviews in Mineralogy and Geochemistry", volume="83", pages="535-575",
       doi="10.2138/rmg.2017.83.12"),
    _r(key="chakraborty2010", authors="Chakraborty, S.", year=2010,
       title="Diffusion coefficients in olivine, wadsleyite and ringwoodite",
       journal="Reviews in Mineralogy and Geochemistry", volume="72", pages="603-639",
       doi="10.2138/rmg.2010.72.13"),
    _r(key="cherniak2010", authors="Cherniak, D. J.", year=2010,
       title="Cation diffusion in feldspars",
       journal="Reviews in Mineralogy and Geochemistry", volume="72", pages="691-733",
       doi="10.2138/rmg.2010.72.15"),
    _r(key="cherniak_dimanov2010", authors="Cherniak, D. J. and Dimanov, A.", year=2010,
       title="Diffusion in pyroxene, mica and amphibole",
       journal="Reviews in Mineralogy and Geochemistry", volume="72", pages="641-690",
       doi="10.2138/rmg.2010.72.14"),
    _r(key="vanorman_crispin2010", authors="Van Orman, J. A. and Crispin, K. L.", year=2010,
       title="Diffusion in oxides",
       journal="Reviews in Mineralogy and Geochemistry", volume="72", pages="757-825",
       doi="10.2138/rmg.2010.72.17"),
    _r(key="lasaga1983", authors="Lasaga, A. C.", year=1983,
       title="Geospeedometry: an extension of geothermometry",
       journal="In: Saxena, S. K. (ed) Kinetics and Equilibrium in Mineral Reactions. Springer, New York",
       pages="81-114", kind="incollection", doi="10.1007/978-1-4612-5587-1_3"),
    _r(key="frost1991", authors="Frost, B. R.", year=1991,
       title="Introduction to oxygen fugacity and its petrologic importance",
       journal="Reviews in Mineralogy", volume="25", pages="1-9"),
    _r(key="deer1992", authors="Deer, W. A. and Howie, R. A. and Zussman, J.", year=1992,
       title="An Introduction to the Rock-Forming Minerals, 2nd edition",
       journal="Longman, Harlow", pages="696 pp", kind="book"),
    _r(key="dmg2025", authors="Dohmen, R. and Chakraborty, S. and course lecturers", year=2025,
       title="DMG Short Course: Diffusion in minerals -- lectures, practicals and supplementary scripts (Diffusion_equation_for_plag.pdf, Diff_Model_Sr_in_Plag_implicit.m, MCdiff_OlFo_MO.m)",
       journal="Ruhr-Universitaet Bochum, October 2025", kind="misc",
       note="unpublished course material, local copy in Volcanology/Diffusion Chronometry/presentations"),
    # ---- physical constants ---------------------------------------------------
    _r(key="codata2018", authors="Tiesinga, E. and Mohr, P. J. and Newell, D. B. and Taylor, B. N.", year=2021,
       title="CODATA recommended values of the fundamental physical constants: 2018",
       journal="Reviews of Modern Physics", volume="93", pages="025010", doi="10.1103/RevModPhys.93.025010"),
    _r(key="iupac2021", authors="Prohaska, T. and Irrgeher, J. and Benefield, J. and others", year=2022,
       title="Standard atomic weights of the elements 2021 (IUPAC Technical Report)",
       journal="Pure and Applied Chemistry", volume="94", pages="573-600", doi="10.1515/pac-2019-0603"),
    _r(key="iau_julian_year", authors="International Astronomical Union", year=1976,
       title="Resolution on the Julian year of 365.25 days (IAU General Assembly XVI, Grenoble)", kind="misc"),
    # ---- oxygen buffers -------------------------------------------------------
    _r(key="oneill_pownceby1993", authors="O'Neill, H. St. C. and Pownceby, M. I.", year=1993,
       title="Thermodynamic data from redox reactions at high temperatures. I. An experimental and theoretical assessment of the electrochemical method using stabilized zirconia electrolytes, with revised values for the Fe-'FeO', Co-CoO, Ni-NiO and Cu-Cu2O oxygen buffers, and new data for the W-WO2 buffer",
       journal="Contributions to Mineralogy and Petrology", volume="114", pages="296-314", doi="10.1007/BF01046533"),
    _r(key="oneill1987", authors="O'Neill, H. St. C.", year=1987,
       title="Quartz-fayalite-iron and quartz-fayalite-magnetite equilibria and the free energy of formation of fayalite (Fe2SiO4) and magnetite (Fe3O4)",
       journal="American Mineralogist", volume="72", pages="67-75"),
    # ---- diffusion modelling methodology ------------------------------------
    _r(key="costa2003", authors="Costa, F. and Chakraborty, S. and Dohmen, R.", year=2003,
       title="Diffusion coupling between trace and major elements and a model for calculation of magma residence times using plagioclase",
       journal="Geochimica et Cosmochimica Acta", volume="67", pages="2189-2200", doi="10.1016/S0016-7037(02)01345-5"),
    _r(key="costa_chakraborty2004", authors="Costa, F. and Chakraborty, S.", year=2004,
       title="Decadal time gaps between mafic intrusion and silicic eruption obtained from chemical zoning patterns in olivine",
       journal="Earth and Planetary Science Letters", volume="227", pages="517-530", doi="10.1016/j.epsl.2004.08.011"),
    _r(key="dohmen_blundy2014", authors="Dohmen, R. and Blundy, J.", year=2014,
       title="A predictive thermodynamic model for element partitioning between plagioclase and melt as a function of pressure, temperature and composition",
       journal="American Journal of Science", volume="314", pages="1319-1372", doi="10.2475/09.2014.04"),
    _r(key="bindeman1998", authors="Bindeman, I. N. and Davis, A. M. and Drake, M. J.", year=1998,
       title="Ion microprobe study of plagioclase-basalt partition experiments at natural concentration levels of trace elements",
       journal="Geochimica et Cosmochimica Acta", volume="62", pages="1175-1193", doi="10.1016/S0016-7037(98)00047-7"),
    _r(key="zellmer1999", authors="Zellmer, G. F. and Blake, S. and Vance, D. and Hawkesworth, C. and Turner, S.", year=1999,
       title="Plagioclase residence times at two island arc volcanoes (Kameni Islands, Santorini, and Soufriere, St. Vincent) determined by Sr diffusion systematics",
       journal="Contributions to Mineralogy and Petrology", volume="136", pages="345-357", doi="10.1007/s004100050543"),
    _r(key="petrone2016", authors="Petrone, C. M. and Bugatti, G. and Braschi, E. and Tommasini, S.", year=2016,
       title="Pre-eruptive magmatic processes re-timed using a non-isothermal approach to magma chamber dynamics",
       journal="Nature Communications", volume="7", pages="12946", doi="10.1038/ncomms12946"),
    _r(key="saunders2012", authors="Saunders, K. and Blundy, J. and Dohmen, R. and Cashman, K.", year=2012,
       title="Linking petrology and seismology at an active volcano",
       journal="Science", volume="336", pages="1023-1027", doi="10.1126/science.1220066"),
    _r(key="tomiya2013", authors="Tomiya, A. and Miyagi, I. and Saito, G. and Geshi, N.", year=2013,
       title="Short time scales of magma-mixing processes prior to the 2011 eruption of Shinmoedake volcano, Kirishima volcanic group, Japan",
       journal="Bulletin of Volcanology", volume="75", pages="750", doi="10.1007/s00445-013-0750-1"),
    _r(key="bradshaw_kent2017", authors="Bradshaw, R. W. and Kent, A. J. R.", year=2017,
       title="The analytical limits of modeling short diffusion timescales",
       journal="Chemical Geology", volume="466", pages="667-677", doi="10.1016/j.chemgeo.2017.07.018"),
    _r(key="shea2015", authors="Shea, T. and Lynn, K. J. and Garcia, M. O.", year=2015,
       title="Cracking the olivine zoning code: Distinguishing between crystal growth and diffusion",
       journal="Geology", volume="43", pages="935-938", doi="10.1130/G37082.1"),
    _r(key="krimer_costa2017", authors="Krimer, D. and Costa, F.", year=2017,
       title="Evaluation of the effects of 3D diffusion, crystal geometry, and initial conditions on retrieved time-scales from Fe-Mg zoning in natural oriented orthopyroxene crystals",
       journal="Geochimica et Cosmochimica Acta", volume="196", pages="271-288", doi="10.1016/j.gca.2016.09.037"),
    _r(key="ganguly1988", authors="Ganguly, J. and Bhattacharya, R. N. and Chakraborty, S.", year=1988,
       title="Convolution effect in the determination of compositional profiles and diffusion coefficients by microprobe step scans",
       journal="American Mineralogist", volume="73", pages="901-909"),
    # ---- olivine ---------------------------------------------------------------
    _r(key="dohmen_chakraborty2007", authors="Dohmen, R. and Chakraborty, S.", year=2007,
       title="Fe-Mg diffusion in olivine II: point defect chemistry, change of diffusion mechanisms and a model for calculation of diffusion coefficients in natural olivine",
       journal="Physics and Chemistry of Minerals", volume="34", pages="409-430", doi="10.1007/s00269-007-0158-6",
       note="see also erratum: Phys Chem Minerals 34:597-598, doi:10.1007/s00269-007-0185-3"),
    _r(key="dohmen2007", authors="Dohmen, R. and Becker, H.-W. and Chakraborty, S.", year=2007,
       title="Fe-Mg diffusion in olivine I: experimental determination between 700 and 1,200 C as a function of composition, crystal orientation and oxygen fugacity",
       journal="Physics and Chemistry of Minerals", volume="34", pages="389-407", doi="10.1007/s00269-007-0157-7"),
    _r(key="chakraborty1997", authors="Chakraborty, S.", year=1997,
       title="Rates and mechanisms of Fe-Mg interdiffusion in olivine at 980-1300 C",
       journal="Journal of Geophysical Research", volume="102", pages="12317-12331", doi="10.1029/97JB00208"),
    _r(key="petry2004", authors="Petry, C. and Chakraborty, S. and Palme, H.", year=2004,
       title="Experimental determination of Ni diffusion coefficients in olivine and their dependence on temperature, composition, oxygen fugacity, and crystallographic orientation",
       journal="Geochimica et Cosmochimica Acta", volume="68", pages="4179-4188", doi="10.1016/j.gca.2004.02.024"),
    _r(key="holzapfel2007", authors="Holzapfel, C. and Chakraborty, S. and Rubie, D. C. and Frost, D. J.", year=2007,
       title="Effect of pressure on Fe-Mg, Ni and Mn diffusion in (FexMg1-x)2SiO4 olivine",
       journal="Physics of the Earth and Planetary Interiors", volume="162", pages="186-198", doi="10.1016/j.pepi.2007.04.009"),
    _r(key="coogan2005ca", authors="Coogan, L. A. and Hain, A. and Stahl, S. and Chakraborty, S.", year=2005,
       title="Experimental determination of the diffusion coefficient for calcium in olivine between 900 C and 1500 C",
       journal="Geochimica et Cosmochimica Acta", volume="69", pages="3683-3694", doi="10.1016/j.gca.2005.03.002"),
    # ---- orthopyroxene ---------------------------------------------------------
    _r(key="dohmen2016", authors="Dohmen, R. and Ter Heege, J. H. and Becker, H.-W. and Chakraborty, S.", year=2016,
       title="Fe-Mg interdiffusion in orthopyroxene",
       journal="American Mineralogist", volume="101", pages="2210-2221", doi="10.2138/am-2016-5815"),
    _r(key="ganguly_tazzoli1994", authors="Ganguly, J. and Tazzoli, V.", year=1994,
       title="Fe2+-Mg interdiffusion in orthopyroxene: Retrieval from the data on intracrystalline exchange reaction",
       journal="American Mineralogist", volume="79", pages="930-937"),
    _r(key="schwandt1998", authors="Schwandt, C. S. and Cygan, R. T. and Westrich, H. R.", year=1998,
       title="Magnesium self-diffusion in orthoenstatite",
       journal="Contributions to Mineralogy and Petrology", volume="130", pages="390-396", doi="10.1007/s004100050373"),
    # ---- clinopyroxene ---------------------------------------------------------
    _r(key="muller2013", authors="Mueller, T. and Dohmen, R. and Becker, H. W. and ter Heege, J. H. and Chakraborty, S.", year=2013,
       title="Fe-Mg interdiffusion rates in clinopyroxene: experimental data and implications for Fe-Mg exchange geothermometers",
       journal="Contributions to Mineralogy and Petrology", volume="166", pages="1563-1576", doi="10.1007/s00410-013-0941-y"),
    _r(key="dimanov_wiedenbeck2006", authors="Dimanov, A. and Wiedenbeck, M.", year=2006,
       title="(Fe,Mn)-Mg interdiffusion in natural diopside: effect of pO2",
       journal="European Journal of Mineralogy", volume="18", pages="705-718", doi="10.1127/0935-1221/2006/0018-0705"),
    _r(key="dimanov_sautter2000", authors="Dimanov, A. and Sautter, V.", year=2000,
       title="'Average' interdiffusion of (Fe,Mn)-Mg in natural diopside",
       journal="European Journal of Mineralogy", volume="12", pages="749-760", doi="10.1127/ejm/12/4/0749"),
    _r(key="brady_mccallister1983", authors="Brady, J. B. and McCallister, R. H.", year=1983,
       title="Diffusion data for clinopyroxenes from homogenization and self-diffusion experiments",
       journal="American Mineralogist", volume="68", pages="95-105"),
    _r(key="vanorman2001", authors="Van Orman, J. A. and Grove, T. L. and Shimizu, N.", year=2001,
       title="Rare earth element diffusion in diopside: influence of temperature, pressure, and ionic radius, and an elastic model for diffusion in silicates",
       journal="Contributions to Mineralogy and Petrology", volume="141", pages="687-703", doi="10.1007/s004100100269"),
    _r(key="sneeringer1984", authors="Sneeringer, M. and Hart, S. R. and Shimizu, N.", year=1984,
       title="Strontium and samarium diffusion in diopside",
       journal="Geochimica et Cosmochimica Acta", volume="48", pages="1589-1608", doi="10.1016/0016-7037(84)90415-6"),
    _r(key="coogan2005li", authors="Coogan, L. A. and Kasemann, S. A. and Chakraborty, S.", year=2005,
       title="Rates of hydrothermal cooling of new oceanic upper crust derived from lithium-geospeedometry",
       journal="Earth and Planetary Science Letters", volume="240", pages="415-424", doi="10.1016/j.epsl.2005.09.020"),
    # ---- plagioclase -----------------------------------------------------------
    _r(key="vanorman2014", authors="Van Orman, J. A. and Cherniak, D. J. and Kita, N. T.", year=2014,
       title="Magnesium diffusion in plagioclase: Dependence on composition, and implications for thermal resetting of the 26Al-26Mg early solar system chronometer",
       journal="Earth and Planetary Science Letters", volume="385", pages="79-88", doi="10.1016/j.epsl.2013.10.026"),
    _r(key="faak2013", authors="Faak, K. and Chakraborty, S. and Coogan, L. A.", year=2013,
       title="Mg in plagioclase: Experimental calibration of a new geothermometer and diffusion coefficients",
       journal="Geochimica et Cosmochimica Acta", volume="123", pages="195-217", doi="10.1016/j.gca.2013.05.009"),
    _r(key="giletti_casserly1994", authors="Giletti, B. J. and Casserly, J. E. D.", year=1994,
       title="Strontium diffusion kinetics in plagioclase feldspars",
       journal="Geochimica et Cosmochimica Acta", volume="58", pages="3785-3793", doi="10.1016/0016-7037(94)90363-8"),
    _r(key="cherniak_watson1994", authors="Cherniak, D. J. and Watson, E. B.", year=1994,
       title="A study of strontium diffusion in plagioclase using Rutherford backscattering spectroscopy",
       journal="Geochimica et Cosmochimica Acta", volume="58", pages="5179-5190", doi="10.1016/0016-7037(94)90303-4"),
    _r(key="cherniak2002", authors="Cherniak, D. J.", year=2002,
       title="Ba diffusion in feldspar",
       journal="Geochimica et Cosmochimica Acta", volume="66", pages="1641-1650", doi="10.1016/S0016-7037(01)00866-3"),
    _r(key="cherniak1996", authors="Cherniak, D. J.", year=1996,
       title="Strontium diffusion in sanidine and albite, and general comments on strontium diffusion in alkali feldspars",
       journal="Geochimica et Cosmochimica Acta", volume="60", pages="5037-5043",
       doi="10.1016/S0016-7037(96)00293-1"),
    _r(key="cherniak_watson2020", authors="Cherniak, D. J. and Watson, E. B.", year=2020,
       title="Ti diffusion in feldspar",
       journal="American Mineralogist", volume="105", pages="1040-1051",
       doi="10.2138/am-2020-7272", note="open access (CC-BY)"),
    _r(key="giletti_shanahan1997", authors="Giletti, B. J. and Shanahan, T. M.", year=1997,
       title="Alkali diffusion in plagioclase feldspar",
       journal="Chemical Geology", volume="139", pages="3-20", doi="10.1016/S0009-2541(97)00026-0"),
    _r(key="grove1984", authors="Grove, T. L. and Baker, M. B. and Kinzler, R. J.", year=1984,
       title="Coupled CaAl-NaSi diffusion in plagioclase feldspar: Experiments and applications to cooling rate speedometry",
       journal="Geochimica et Cosmochimica Acta", volume="48", pages="2113-2121", doi="10.1016/0016-7037(84)90391-6"),
    _r(key="liu_yund1992", authors="Liu, M. and Yund, R. A.", year=1992,
       title="NaSi-CaAl interdiffusion in plagioclase",
       journal="American Mineralogist", volume="77", pages="275-283"),
    # ---- magnetite -------------------------------------------------------------
    _r(key="freer_hauptman1978", authors="Freer, R. and Hauptman, Z.", year=1978,
       title="An experimental study of magnetite-titanomagnetite interdiffusion",
       journal="Physics of the Earth and Planetary Interiors", volume="16", pages="223-231", doi="10.1016/0031-9201(78)90015-8"),
    _r(key="aragon1984", authors="Aragon, R. and McCallister, R. H. and Harrison, H. R.", year=1984,
       title="Cation diffusion in titanomagnetites",
       journal="Contributions to Mineralogy and Petrology", volume="85", pages="174-185", doi="10.1007/BF00371707"),
    _r(key="sievwright2020", authors="Sievwright, R. H. and O'Neill, H. St. C. and Tolley, J. and Wilkinson, J. J. and Berry, A. J.", year=2020,
       title="Diffusion and partition coefficients of minor and trace elements in magnetite as a function of oxygen fugacity at 1150 C",
       journal="Contributions to Mineralogy and Petrology", volume="175", pages="40", doi="10.1007/s00410-020-01679-z",
       note="local library copy: sievwright_etal_2020_magnetite_trace_diffusion.pdf. Table 5 gives the vacancy and interstitial constants at 1150 C"),
    # ---- applications and secondary sources cited by the registry ----------
    _r(key="sato2022", authors="Sato, E. and Ban, M. and Yoshida, T. and Andrews, B.", year=2022,
       title="Magma plumbing system and eruption processes of the Okama pyroclastics, Zao volcano, revealed by orthopyroxene Fe-Mg diffusion chronometry",
       journal="Journal of Volcanology and Geothermal Research", volume="429", pages="107607",
       doi="10.1016/j.jvolgeores.2022.107607",
       note="local library copy: sato_etal_2022_zao_opx_diffusion.pdf"),
    _r(key="polo_sanchez2023", authors="Polo-Sanchez, A. and Druitt, T. H. and Cluzel, N. and Devidal, J.-L.", year=2023,
       title="Pyroxene diffusion chronometry of the magmatic plumbing system of Santorini volcano",
       journal="Frontiers in Earth Science", volume="11", pages="1149446",
       note="local library copy: polo-sanchez_etal_2023_santorini_pyroxene_diffusion.pdf"),
    _r(key="ostorero2022", authors="Ostorero, L. and Balcone-Boissard, H. and Boudon, G. and others", year=2022,
       title="Correlated petrology and seismicity indicate rapid magma accumulation prior to eruption of Kizimen volcano, Kamchatka",
       journal="Communications Earth & Environment", volume="3", pages="290",
       doi="10.1038/s43247-022-00622-3"),
    _r(key="mutch2021", authors="Mutch, E. J. F. and Maclennan, J. and Madden-Nadeau, A. L.", year=2021,
       title="DFENS: Diffusion chronometry using finite elements and nested sampling",
       journal="Geochemistry, Geophysics, Geosystems", volume="22", pages="e2020GC009303",
       doi="10.1029/2020GC009303",
       note="source code at https://github.com/EuanMutch/DFENS. It samples the diffusion-law parameters from their covariance matrix"),
    _r(key="hartley2016", authors="Hartley, M. E. and Morgan, D. J. and Maclennan, J. and Edmonds, M. and Thordarson, T.", year=2016,
       title="Tracking timescales of short-term precursors to large basaltic fissure eruptions through Fe-Mg diffusion in olivine",
       journal="Earth and Planetary Science Letters", volume="439", pages="58-70",
       doi="10.1016/j.epsl.2016.01.018"),
    _r(key="girona_costa2013", authors="Girona, T. and Costa, F.", year=2013,
       title="DIPRA: A user-friendly program to model multi-element diffusion in olivine with applications to timescales of magmatic processes",
       journal="Geochemistry, Geophysics, Geosystems", volume="14", pages="422-431",
       doi="10.1029/2012GC004427"),
    _r(key="aggarwal_dieckmann2002", authors="Aggarwal, S. and Dieckmann, R.", year=2002,
       title="Point defects and cation tracer diffusion in (Ti_x Fe_(1-x))_(3-d) O_4. II. Cation tracer diffusion",
       journal="Physics and Chemistry of Minerals", volume="29", pages="707-718",
       doi="10.1007/s00269-002-0282-2",
       note="primary source of the magnetite tracer data tabulated by Van Orman & Crispin (2010)"),
    _r(key="dieckmann1987", authors="Dieckmann, R. and Mason, T. O. and Hodge, J. D. and Schmalzried, H.", year=1978,
       title="Defects and cation diffusion in magnetite (III): tracer diffusion of foreign tracer cations as a function of temperature and oxygen potential",
       journal="Berichte der Bunsengesellschaft fuer physikalische Chemie", volume="82", pages="778-783",
       note="primary source of the Cr and Al entries in Van Orman & Crispin (2010) Table 12"),
    _r(key="huebner1971", authors="Huebner, J. S.", year=1971,
       title="Buffering techniques for hydrostatic systems at elevated pressures",
       journal="In: Ulmer, G. C. (ed) Research Techniques for High Pressure and High Temperature. Springer",
       pages="123-177", kind="incollection",
       note="buffer equations used by Van Orman & Crispin (2010) to compute their Tables 10 and 11"),
    _r(key="allan2013", authors="Allan, A. S. R. and Morgan, D. J. and Wilson, C. J. N. and Millet, M.-A.", year=2013,
       title="From mush to eruption in centuries: assembly of the super-sized Oruanui magma body",
       journal="Contributions to Mineralogy and Petrology", volume="166", pages="143-164",
       doi="10.1007/s00410-013-0869-2",
       note="source of the fO2-dependent form of the Ganguly & Tazzoli (1994) orthopyroxene law"),
    _r(key="latourrette_wasserburg1998", authors="LaTourrette, T. and Wasserburg, G. J.", year=1998,
       title="Mg diffusion in anorthite: implications for the formation of early solar system planetesimals",
       journal="Earth and Planetary Science Letters", volume="158", pages="91-108",
       doi="10.1016/S0012-821X(98)00048-X"),
    _r(key="dias2025", authors="Dias, M. A. and Dohmen, R. and Behrens, H.", year=2025,
       title="Fe-Mg interdiffusion in orthopyroxene: complex interdependencies of temperature, composition and oxygen fugacity",
       journal="Geochimica et Cosmochimica Acta", volume="395", pages="195-211",
       doi="10.1016/j.gca.2025.03.002",
       note="local library copy: dias_dohmen_2025_opx_femg_interdependencies.pdf. Eqs 22-25 give separate parameterisations above and below log fO2 = -10 Pa"),
    _r(key="dias_dohmen2024", authors="Dias, M. A. and Dohmen, R.", year=2024,
       title="Experimental determination of Fe-Mg interdiffusion in orthopyroxene as a function of Fe content",
       journal="Contributions to Mineralogy and Petrology", volume="179", pages="36",
       doi="10.1007/s00410-024-02110-7",
       note="local library copy: dias_dohmen_2024_opx_fe_content_diffusion.pdf"),
    _r(key="dias2025ree", authors="Dias, M. A. and Dohmen, R. and Hartmann, N.", year=2025,
       title="Diffusion of Eu, Ce and Lu in orthopyroxene",
       journal="Geochimica et Cosmochimica Acta", volume="410", pages="85-100",
       doi="10.1016/j.gca.2025.09.038",
       note="local library copy: dias_dohmen_2025_opx_eu_ce_lu.pdf"),
    _r(key="pohl2024", authors="Pohl, F. and Behrens, H. and Oeser, M. and Marxer, F. and Dohmen, R.", year=2024,
       title="Li diffusion in plagioclase crystals and glasses - implications for timescales of geological processes",
       journal="European Journal of Mineralogy", volume="36", pages="985-1003",
       doi="10.5194/ejm-36-985-2024",
       note="open access, local library copy: pohl_etal_2024_plagioclase_li.pdf"),
    _r(key="audetat2026", authors="Audetat, A. and Grocolas, T. and Mutch, E. J. F.", year=2026,
       title="Ti-in-quartz and Sr-Ba-Mg-in-feldspars diffusion chronometry: a review of available diffusion data, and a critical evaluation of applications to natural samples",
       journal="Journal of Petrology", volume="", pages="egag078",
       doi="10.1093/petrology/egag078",
       note="accepted manuscript (advance article), page numbers not yet assigned. Local library copy: audedat_grocolas_mutch_2025_manuscript_ti_quartz_sr_ba_mg_feldspars.pdf. Eq. 1 parameterises Mg diffusion in plagioclase from Faak et al. (2013) and Van Orman et al. (2014)"),
    _r(key="morgan2004", authors="Morgan, D. J. and Blake, S. and Rogers, N. W. and others", year=2004,
       title="Time scales of crystal residence and magma chamber volume from modelling of diffusion profiles in phenocrysts: Vesuvius 1944",
       journal="Earth and Planetary Science Letters", volume="222", pages="933-946",
       doi="10.1016/j.epsl.2004.03.030"),
    _r(key="kress_carmichael1991", authors="Kress, V. C. and Carmichael, I. S. E.", year=1991,
       title="The compressibility of silicate liquids containing Fe2O3 and the effect of composition, temperature, oxygen fugacity and pressure on their redox states",
       journal="Contributions to Mineralogy and Petrology", volume="108", pages="82-92",
       doi="10.1007/BF00307328"),
    _r(key="grocolas2025", authors="Grocolas, T. and Bloch, E. M. and Bouvier, A.-S. and Muentener, O.", year=2025,
       title="Diffusion of Sr and Ba in plagioclase: composition and silica activity dependencies, and application to volcanic rocks",
       journal="Earth and Planetary Science Letters", volume="651", pages="119141",
       doi="10.1016/j.epsl.2024.119141",
       note="open access (CC-BY), local library copy: grocolas_etal_2025_plagioclase_sr_ba.pdf. Finds Sr diffusion in plagioclase 1.5-2 orders of magnitude SLOWER than Giletti & Casserly (1994). Ba diffusion is similar to earlier work. Experiments on oligoclase and labradorite, 900-1200 C, 1 atm, with aSiO2 buffered, and no resolvable dependence on aSiO2 or crystal orientation. Eqs 7-8 are implemented, and eqs 12-14 (their Monte Carlo re-fits of the older data) are used for the older entries"),
    _r(key="grocolas2025cmp", authors="Grocolas, T. and Muentener, O. and Bloch, E. M. and Escrig, S. and Ulyanov, A. and Bouvier, A.-S.", year=2025,
       title="Cooling rates and melt extraction timescales determined by diffusion chronometry on shallow crustal plutonic rocks",
       journal="Contributions to Mineralogy and Petrology", volume="180", pages="45",
       doi="10.1007/s00410-025-02238-0",
       note="open access, PMC12254178. Applies the Grocolas et al. (2025) Sr and Ba diffusivities to the Adamello batholith"),
    _r(key="chamberlain2014", authors="Chamberlain, K. J. and Morgan, D. J. and Wilson, C. J. N.", year=2014,
       title="Timescales of mixing and mobilisation in the Bishop Tuff magma body: perspectives from diffusion chronometry",
       journal="Contributions to Mineralogy and Petrology", volume="168", pages="1034",
       doi="10.1007/s00410-014-1034-2",
       note="local library copy: chamberlain_etal_2014_bishop_tuff_diffusion.pdf. Feldspar microprobe traverses with a 5 um defocused beam"),
    _r(key="druitt2012", authors="Druitt, T. H. and Costa, F. and Deloule, E. and Dungan, M. and Scaillet, B.", year=2012,
       title="Decadal to monthly timescales of magma transfer and reservoir growth at a caldera volcano",
       journal="Nature", volume="482", pages="77-80", doi="10.1038/nature10706",
       note="Supplementary Table 1 gives measured An, Mg, Sr, Ba, Ti, K, La and Ce profiles across Minoan plagioclase phenocrysts. Diffusor ships crystal S82-30A 12 as a worked example"),
]}


def get(key: str) -> Reference:
    try:
        return REFERENCES[key]
    except KeyError as exc:
        raise KeyError(f"Unknown citation key '{key}'. Add it to diffusor/references.py.") from exc


def cite(key: str, where: str = "") -> str:
    """Compact citation, e.g. ``cite('crank1975', 'eq. 2.14')`` -> 'Crank (1975), eq. 2.14'."""
    s = get(key).short()
    return f"{s}, {where}" if where else s


def format_reference(key: str) -> str:
    return get(key).full()


def to_bibtex(ref: Reference) -> str:
    fields = {
        "author": ref.authors,
        "year": str(ref.year),
        "title": ref.title,
    }
    if ref.kind == "article":
        fields["journal"] = ref.journal
        if ref.volume:
            fields["volume"] = ref.volume
        if ref.pages:
            fields["pages"] = ref.pages
    elif ref.kind == "book":
        fields["publisher"] = ref.journal
    elif ref.kind == "incollection":
        fields["booktitle"] = ref.journal
        if ref.pages:
            fields["pages"] = ref.pages
    else:
        if ref.journal:
            fields["howpublished"] = ref.journal
    if ref.doi:
        fields["doi"] = ref.doi
    if ref.note:
        fields["note"] = ref.note
    body = ",\n".join(f"  {k} = {{{v}}}" for k, v in fields.items())
    return f"@{ref.kind}{{{ref.key},\n{body}\n}}"


def all_keys():
    return sorted(REFERENCES)
