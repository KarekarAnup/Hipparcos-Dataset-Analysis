# Hipparcos Dataset: Stellar Data Analysis Using Machine Learning

### About Dataset: - 
(https://heasarc.gsfc.nasa.gov/W3Browse/all/hipparcos.html)

The Hipparcos and Tycho Catalogues are the primary products of the European Space Agency's astrometric mission, Hipparcos. The satellite, which operated for four years, returned high quality scientific data from November 1989 to March 1993.
Each of the catalogues contains a large quantity of very high quality astrometric and photometric data. In addition there are associated annexes featuring variability and double/multiple star data, and solar system astrometric and photometric measurements. In the case of the Hipparcos Catalogue, the principal parts are provided in both printed and machine-readable form (on CDROM). In the case of the Tycho Catalogue, results are provided in machine-readable form only (on CDROM). Although in general only the final reduced and calibrated astrometric and photometric data are provided, some auxiliary files containing results from intermediate stages of the data processing, of relevance for the more-specialised user, have also been retained for publication. (Some, but not all, data files are available from the Centre de Donnees astronomiques de Strasbourg.)

The global data analysis tasks, proceeding from nearly 1000 Gbit of raw satellite data to the final catalogues, was a lengthy and complex process, and was undertaken by the NDAC and FAST Consortia, together responsible for the production of the Hipparcos Catalogue, and the Tycho Consortium, responsible for the production of the Tycho Catalogue. A fourth scientific consortium, the INCA Consortium, was responsible for the construction of the Hipparcos observing programme, compiling the best-available data for the selected stars before launch into the Hipparcos Input Catalogue. The production of the Hipparcos and Tycho Catalogues marks the formal end of the involvement in the mission by the European Space Agency and the four scientific consortia.

For much more information about this catalog, such as fuller descriptions of the parameters, the user is urged to check the Hipparcos and Tycho Catalogs website at https://www.cosmos.esa.int/web/hipparcos/catalogues.

### Citation
- The data has been collected by the European Space Agency's astrometric satellite , Hipparcos , active from August 1989 to August 1993 .
- The catalog was published on July 1997 .

### Reference: - 

Perryman M.A.C., Lindegren L., Kovalevsky J., Hog E., Bastian U.,
   Bernacca P.L., Creze M., Donati F., Grenon M., Grewing M.,
   van Leeuwen F., van der Marel H., Mignard F., Murray C.A.,
   Le Poole R.S., Schrijver H., Turon C., Arenou F., Froeschle M.,
   Petersen C.S., "The Hipparcos Catalogue"
   =1997A&A...323L..49P

Lindegren L., Mignard F., Soederhjelm S., Badiali M., Bernstein H.H.,
   Lampens P., Pannunzio R., Arenou F., Bernacca P.L., Falin J.L.,
   Froeschle M., Kovalevsky J., Martin C., Perryman M.A.C., Wielen R.
   "Double star data in the Hipparcos Catalogue"
   =1997A&A...323L..53L

Hog E., Baessgen G., Bastian U., Egret D., Fabricius C., Grossmann V.,
   Halbwachs J.L., Makarov V.V., Perryman M.A.C., Schwekendiek P.,
   Wagner K., Wicenec A., "The Tycho Catalogue"
   =1997A&A...323L..57H

van Leeuwen F., Evans D.W., Grenon M., Grossmann V., Mignard F.,
   Perryman M.A.C., "The Hipparcos mission: photometric data."
   =1997A&A...323L..61V

### Project Description:- 
#### Hipparcos Dataset: Stellar Data Analysis Using Machine Learning

This project involves comprehensive exploration and analysis of a dataset containing astronomical parameters for stars, with a focus on deriving meaningful insights and classifications using various machine learning techniques. The dataset includes key parameters such as visual magnitude (Vmag), parallax (Plx), right ascension (RAdeg), declination (DEdeg), proper motion in right ascension (pmRA), proper motion in declination (pmDE), B-V color index (B-V), and spectral type (SpType). The following key analyses and machine learning applications were performed:

### 1. Exploratory Data Analysis (EDA)
- **Distance Calculation**: The distance to each star was computed using the parallax angle with the formula \(d = \frac{1}{\text{Plx}}\), after converting Plx from milliarcseconds to arcseconds.
- **Absolute Magnitude**: Calculated using the formula \(M = Vmag + 5 \log_{10}(\text{Plx}) + 5\), providing a measure of the star's intrinsic brightness.
- **Proper Motion**: Combined proper motion in RA and DEC to calculate the total proper motion \(\mu = \sqrt{(\text{pmRA} \cdot \cos(\text{DEdeg}))^2 + (\text{pmDE})^2}\).
- **Tangential Velocity**: Computed using the formula \(v_t = 4.74 \times \frac{\mu}{\text{Plx}}\), where \(v_t\) is in km/s.
- **Hertzsprung-Russell Diagram**: Plotted using calculated absolute magnitudes and B-V color indices to visualize stellar properties.
- **Color Index and Temperature**: Analyzed the B-V color index to estimate the surface temperature of stars.
- **Spectral Classification**: Examined spectral types to categorize stars into main sequence, giants, dwarfs, etc.

### 2. Machine Learning Applications
#### Clustering Algorithms
- **K-Means Clustering**:
  - **Application**: Grouped stars with similar properties such as distance, proper motion, and color index to identify potential star clusters.
  - **Implementation**: Used `sklearn.cluster.KMeans` to perform clustering and visualized clusters in a 3D space of RA, DEC, and distance.
- **DBSCAN (Density-Based Spatial Clustering of Applications with Noise)**:
  - **Application**: Identified star associations or clusters in regions with high stellar density, especially effective in noisy datasets.
  - **Implementation**: Applied `sklearn.cluster.DBSCAN` to detect non-linear star clusters.

#### Classification Algorithms
- **Random Forest Classifier**:
  - **Application**: Classified stars into spectral classes using features like B-V color index, visual magnitude, and proper motion.
  - **Implementation**: Utilized `sklearn.ensemble.RandomForestClassifier` to build and train the classifier, evaluating performance with metrics such as accuracy, precision, and recall.
- **Support Vector Machines (SVM)**:
  - **Application**: Predicted the spectral type of stars, effective for both binary and multiclass classification tasks.
  - **Implementation**: Employed `sklearn.svm.SVC` to classify stars, optimizing with different kernels like linear, RBF, and polynomial.

#### Regression Algorithms
- **Linear Regression**:
  - **Application**: Predicted continuous variables such as distance based on visual magnitude and parallax.
  - **Implementation**: Trained a linear regression model using `sklearn.linear_model.LinearRegression`, evaluated with RMSE and R² metrics.

### Conclusion
The project successfully utilized a combination of exploratory data analysis and machine learning techniques to analyze stellar data. By calculating key astronomical parameters, clustering stars based on their properties, and classifying them into spectral types, the project provided deep insights into the characteristics and distributions of stars. These methods can be further extended to more complex datasets and refined with advanced machine learning techniques to enhance astronomical research and discoveries.

### HEASARC Parameter Names
The following is a correlation table of the Hipparcos Catalog fields, the parameter names as implemented by the CDS, and the HEASARC names for these parameters:
Hipparcos   CDS Name    HEASARC Name       Description
Cat. Field

--          * New *    Name               /Catalog Designation
H0          Catalog    * Not Displayed *  /Catalogue (H=Hipparcos)

H1          HIP        HIP_Number         /Identifier (HIP number)

H2          Proxy      Prox_10asec        /Proximity flag

H3          RAhms      RA                 /RA in h m s, ICRS (J1991.25)

H4          DEdms      Dec                /Dec in deg ' ", ICRS (J1991.25)

H5          Vmag       Vmag               /Magnitude in Johnson V

H6          VarFlag    Var_Flag           /Coarse variability flag

H7          r_Vmag     Vmag_Source        /Source of magnitude

H8          RAdeg      RA_Deg             /RA in degrees (ICRS, Epoch-J1991.25)

H9          DEdeg      Dec_Deg            /Dec in degrees (ICRS, Epoch-J1991.25)

H10         AstroRef   Astrom_Ref_Dbl     /Reference flag for astrometry

H11         Plx        Parallax           /Trigonometric parallax

H12         pmRA       pm_RA              /Proper motion in RA

H13         pmDE       pm_Dec             /Proper motion in Dec

H14         e_RAdeg    RA_Error           /Standard error in RA*cos(Dec_Deg)

H15         e_DEdeg    Dec_Error          /Standard error in Dec_Deg

H16         e_Plx      Parallax_Error     /Standard error in Parallax

H17         e_pmRA     pm_RA_Error        /Standard error in pmRA

H18         e_pmDE     pm_Dec_Error       /Standard error in pmDE

H19         DE:RA      Crl_Dec_RA         /(DE over RA)xCos(delta)

H20         Plx:RA     Crl_Plx_RA         /(Plx over RA)xCos(delta)

H21         Plx:DE     Crl_Plx_Dec        /(Plx over DE)

H22         pmRA:RA    Crl_pmRA_RA        /(pmRA over RA)xCos(delta)

H23         pmRA:DE    Crl_pmRA_Dec       /(pmRA over DE)


H24         pmRA:Plx   Crl_pmRA_Plx       /(pmRA over Plx)

H25         pmDE:RA    Crl_pmDec_RA       /(pmDE over RA)xCos(delta)

H26         pmDE:DE    Crl_pmDec_Dec      /(pmDE over DE)

H27         pmDE:Plx   Crl_pmDec_Plx      /(pmDE over Plx)

H28         pmDE:pmRA  Crl_pmDec_pmRA     /(pmDE over pmRA)

H29         F1         Reject_Percent     /Percentage of rejected data

H30         F2         Quality_Fit        /Goodness-of-fit parameter

H31         ---        * Not Displayed *  /HIP number (repetition)

H32         BTmag      BT_Mag             /Mean BT magnitude

H33         e_BTmag    BT_Mag_Error       /Standard error on BTmag

H34         VTmag      VT_Mag             /Mean VT magnitude

H35         e_VTmag    VT_Mag_Error       /Standard error on VTmag

H36         m_BTmag    BT_Mag_Ref_Dbl     /Reference flag for BT and VTmag

H37         B-V        BV_Color           /Johnson BV colour

H38         e_B-V      BV_Color_Error     /Standard error on BV

H39         r_B-V      BV_Mag_Source      /Source of BV from Ground or Tycho

H40         V-I        VI_Color           /Colour index in Cousins' system

H41         e_V-I      VI_Color_Error     /Standard error on VI

H42         r_V-I      VI_Color_Source    /Source of VI

H43         CombMag    Mag_Ref_Dbl        /Flag for combined Vmag, BV, VI

H44         Hpmag      Hip_Mag            /Median magnitude in Hipparcos system

H45         e_Hpmag    Hip_Mag_Error      /Standard error on Hpmag

H46         Hpscat     Scat_Hip_Mag       /Scatter of Hpmag

H47         o_Hpmag    N_Obs_Hip_Mag      /Number of observations for Hpmag

H48         m_Hpmag    Hip_Mag_Ref_Dbl    /Reference flag for Hpmag

H49         Hpmax      Hip_Mag_Max        /Hpmag at maximum (5th percentile)

H50         HPmin      Hip_Mag_Min        /Hpmag at minimum (95th percentile)

H51         Period     Var_Period         /Variability period (days)

H52         HvarType   Hip_Var_Type       /Variability type

H53         moreVar    Var_Data_Annex     /Additional data about variability

H54         morePhoto  Var_Curv_Annex     /Light curve Annex

H55         CCDM       CCDM_Id            /CCDM identifier

H56         n_CCDM     CCDM_History       /Historical status flag

H57         Nsys       CCDM_N_Entries     /Number of entries with same CCDM

H58         Ncomp      CCDM_N_Comp        /Number of components in this entry

H59         MultFlag   Dbl_Mult_Annex     /Double and or Multiple Systems flag

H60         Source     Astrom_Mult_Source /Astrometric source flag

H61         Qual       Dbl_Soln_Qual      /Solution quality flag

H62         m_HIP      Dbl_Ref_ID         /Component identifiers

H63         theta      Dbl_Theta          /Position angle between components

H64         rho        Dbl_Rho            /Angular separation of components

H65         e_rho      Rho_Error          /Standard error of rho

H66         dHp        Diff_Hip_Mag       /Magnitude difference of components

H67         e_dHp      dHip_Mag_Error     /Standard error in dHp

H68         Survey     Survey_Star        /Flag indicating a Survey Star

H69         Chart      ID_Chart           /Identification Chart

H70         Notes      Notes              /Existence of notes

H71         HD         HD_Id              /HD number <III 135>

H72         BD         BD_Id              /Bonner DM <I 119>, <I 122>

H73         CoD        CoD_Id             /Cordoba Durchmusterung (DM) <I 114>

H74         CPD        CPD_Id             /Cape Photographic DM <I 108>

H75         (V-I)red   VI_Color_Reduct    /VI used for reductions

H76         SpType     Spect_Type         /Spectral type

H77         r_SpType   Spect_Type_Source  /Source of spectral type

--          * New *    Class              /HEASARC BROWSE classification

##### Parameters
Name: 
Name of the star in the recommended format for Hipparcos stars, as created by concatenating the prefix 'HIP ' and the Hip_Number 
identifier in the original catalog. Entries in the Hipparcos (HIP) Catalog have exactly the same identifier as in the Hipparcos 
Input Catalog (HIC), notice.



RA: 
Right ascension in the specified equinox for epoch J1991.25. This was given in the ICRS reference system (J2000 equator) in the original Hipparcos Catalog, and thus equinox 2000 should be specified to avoid inaccuracies due to the non-rigorous HEASARC coordinate precession algorithm. This parameter was given to a truncated precision of 0.01 seconds of time in the original Hipparcos Catalog. If the 'precise' RA is desired, one should use the value of the parameter RA_deg which contains the complete RA in decimal degrees.

Dec: 
Declination in specified equinox for epoch J1991.25. This was given in the ICRS reference system (J2000 equator) in the original Hipparcos Catalog, and thus equinox 2000 should be specified to avoid inaccuracies due to the non-rigorous HEASARC coordinate precession algorithm. This parameter was given to a truncated precision of 0.1 arcseconds in the original Hipparcos Catalog. If the 'precise' declination is desired, one should use the value of the parameter Dec_deg which contains the complete declination in decimal degrees.

LII: 
Galactic longitude.

BII: 
Galactic latitude.

HIP_Number: 
The Hipparcos Catalog running number, which is the same as the that in the Hipparcos Input Catalog. The star entries are, with a few exceptions, ordered by increasing HIP number, which basically follows the order of the object's right ascension (Equinox J2000) independent of declination.

Prox_10asec: 
A proximity flag which provides a coarse indication of the presence of nearby objects within 10 arcseconds of the position of the given star. If non-blank, it indicates that there are one or more distinct Hipparcos ('H') or Tycho ('T') Catalog entries; if both 'H' and 'T' apply, then 'H' is the adopted value, notice.

Vmag: 
The magnitude in Johnson V band, given to a precision of 0.01 magnitudes in the original Hipparcos Catalog.

Var_Flag: 
A coarse variability flag which indicates if the entry (or one of the components in the case of a multiple system) is variable in its Hipparcos magnitude Hip_mag at the level of:

       1: < 0.06mag ; 2: 0.06-0.6mag ; 3: >0.6mag
  
Vmag_Source: 
The source of the V magnitude:

       G:  ground-based multicolor photometry, either directly in or
           reduced to the Johnson UBV system
       H:  Hipparcos magnitude Hip_mag, combined with information on the
           color index (either V-I or BT_mag-VT_mag), in combination with
           the luminosoty class
       T:  Tycho photometry, i.e., VT_mag and BT_mag-VT_mag
        :  no data available
  
RA_Deg: 
The right ascension expressed in degrees for epoch J1991.25 (JD2448349.0625 (TT)) in the ICRS (International Celestial Reference System, consistent with J2000) reference system, and given to a precision of 10-8 degrees in the original Hipparcos Catalog. There are 263 cases where these fields are missing (no astrometric solution could be found).

Dec_Deg: 
The declination expressed in degrees for epoch J1991.25 (JD2448349.0625 (TT)) in the ICRS (International Celestial Reference System, consistent with J2000) reference system, and given to a precision of 10-8 degrees in the original Hipparcos Catalog. There are 263 cases where these fields are missing (no astrometric solution could be found)

Astrom_Ref_Dbl: 
Reference flag for astrometric parameters of double and multiple systems. This flag indicates that the astrometric parameters refer to:

    A, B etc: the letter indicates the specified component of a double
              or multiple system
           *: the photocentre of a double or multiple system included in
              Part C of the Double and Multiple Systems Annex
           +: the centre of mass: for such an entry, an orbit is given in
              Part O of the Double and Multiple Systems Annex
  
Parallax: 
The trigonometric parallax pi in units of milliarcseconds: thus to calculate the distance D in parsecs, D = 1000/pi. The estimated parallax is given for every star, even if it appears to be insignificant or negative.

PM_RA: 
The proper motion component in the RA direction expressed in milliarcseconds per Julian year (mas/yr), and given with respect to the ICRS reference system: mu_RA* = mu_RA x cos (declination).

PM_Dec: 
The proper motion component in the declination direction expressed in milliarcseconds per Julian year (mas/yr), and given with respect to the ICRS reference system.

RA_Error: 
The standard error in the Right Ascension given at the catalog epoch, J1991.25, and expressed in milliarcseconds: sigma_RA* = sigma_RA x cos (declination).

Dec_Error: 
The standard error in the declination given at the catalog epoch, J1991.25, and expressed in milliarcseconds.

Parallax_Error: 
The standard error in the parallax given in milliarcseconds.

PM_RA_Error: 
The standard error in the proper motion component in the RA direction expressed in milliarcseconds per Julian year (mas/yr): sigma_mu_RA* = sigma_mu_RA x cos (declination).

PM_Dec_Error: 
The standard error in the proper motion component in the declination direction expressed in milliarcseconds per Julian year (mas/yr), sigma_mu_declination.

Crl_Dec_RA: 
The correlation coefficient expressed as a real numerical value (in the printed catalog this is expressed in per cent, notice): (declination over RA).

Crl_Plx_RA: 
The correlation coefficient expressed as a real numerical value (in the printed catalog this is expressed in per cent, notice): (parallax over RA).

Crl_Plx_Dec: 
The correlation coefficient expressed as a real numerical value (in the printed catalog this is expressed in per cent, notice): (parallax over declination).

Crl_Pmra_RA: 
The correlation coefficient expressed as a real numerical value (in the printed catalog this is expressed in per cent, notice): (proper motion in RA over RA).

Crl_Pmra_Dec: 
The correlation coefficient expressed as a real numerical value (in the printed catalog this is expressed in per cent, notice): (proper motion in RA over declination).

Crl_Pmra_Plx: 
The correlation coefficient expressed as a real numerical value (in the printed catalog this is expressed in per cent, notice): (proper motion in RA over parallax).

Crl_Pmdec_RA: 
The correlation coefficient expressed as a real numerical value (in the printed catalog this is expressed in per cent, notice): (proper motion in declination over RA).

Crl_Pmdec_Dec: 
The correlation coefficient expressed as a real numerical value (in the printed catalog this is expressed in per cent, notice): (proper motion in declination over declination).

Crl_Pmdec_Plx: 
The correlation coefficient expressed as a real numerical value (in the printed catalog this is expressed in per cent, notice): (proper motion in declination over parallax).

Crl_Pmdec_Pmra: 
The correlation coefficient expressed as a real numerical value (in the printed catalog this is expressed in per cent, notice): (proper motion in declination over proper motion in RA).

Reject_Percent: 
The percentage of data that had to be rejected in order to obtain an acceptable solution.

Quality_Fit: 
The goodness-of-fit statistic: this number indicates the goodness of fit of the astrometric solution to the accepted data (i.e., excluding the rejected data). For good fits, this should approximately follow a normal distribution with zero mean value and unit standard deviation. Values exceeding, say +3, thus indicate a bad fit to the data.

BT_Mag: 
The mean magnitude in the Tycho photometric system, B_T.

BT_Mag_Error: 
The standard error of the B_T magnitude, BT_mag.

VT_Mag: 
The mean magnitude in the Tycho photometric system, V_T.

VT_Mag_Error: 
The standard error of the V_T magnitude, VT_mag.

BT_Mag_Ref_Dbl: 
a reference flag for BT_mag and VT_mag which indicates, for non-single stars, the component measured in Tycho photometry, or indicates that several components have been directly measured together by Tycho, or have had their Tycho data combined. The flag takes the following values:

   A, B, etc. : the Tycho photometry refers to the designated Hipparcos
                Catalog component
            * : the Tycho photometry refers to all components of the
                relevant Hipparcos entry
            - : the Tycho photometry refers to a single-pointing triple or
                quadruple system, for which only a close pair has been
                observed by Tycho, the other components being too faint
                to be detected by Tycho
  
BV_Color: 
The (B-V) color index in, or reduced to, the Johnson UBV system.

BV_Color_Error: 
The standard error of the (B-V) color index, BV_Color.

BV_Mag_Source: 
The source of the (B-V) color index, BV_Color:

            G: indicates that it was taken from ground-based observations
            T: indicates that it was determined from the transformed Tycho
               (B_T-V_T) data
             : indicates that no data are available
  
VI_Color: 
the (V-I) color index in Cousins' photometric system; it represents the best available (V-I) value at the time of the Hipparcos Catalog publication.

VI_Color_Error: 
The standard error in the (V-I) color index, VI_Color.

VI_Color_Source: 
The Source of the (V-I) color index, VI_Color (see Section 1.3, Appendix 5 of the published Hipparcos Catalog for full details):

       'A'        :for an observation of V-I in Cousins' system;
       'B' to 'K' :when V-I derived from measurements in other
                   bands/photoelectric systems
       'L' to 'P' :when V-I derived from Hipparcos and Star Mapper
                   photometry
       'Q'        :for long-period variables
       'R' to 'T' :when colours are unknown
  
Mag_Ref_Dbl: 
A reference flag for the (B-V) and (V-I) color indices and the V magnitude Vmag (and all their standard errors) which is set to '*' when they refer to the combined light of double or multiple systems which are otherwise resolved by the main mission astrometry and photometry.

Hip_Mag: 
The median magnitude H_P in the Hipparcos photometric system, and defined on the basis of the accepted observations (or field transits) for a given star. Note that the Hipparcos magnitude could not be determined for 14 stars.

Hip_Mag_Error: 
The standard error of the median magnitude H_P.

Scat_Hip_Mag: 
The scatter of the H_P observations.

N_Obs_Hip_Mag: 
The number of H_P observations: this is the number of photometric observations (or field transits) used for the construction of the median, standard error, and scatter in H_P.

Hip_Mag_Ref_Dbl: 
A reference flag for the Hipparcos photometric parameters. For a double or multiple entry, this flag indicates that the photometry refers to:

   A, B, etc. : the specified component of a double or multiple system
            * : combined photometry of a double system, corrected for
                attenuation by the detector's instantaneous field of view
                profile response
            - : combined photometry of a double system, NOT corrected for
                attenuation by the detector's instantaneous field of view
                profile response
  
Hip_Mag_Max: 
The observed magnitude at maximum luminosity. This is defined as the 5th percentile of the epoch photometry.

Hip_Mag_Min: 
The observed magnitude at minimum luminosity. This is defined as the 95th percentile of the epoch photometry.

Var_Period: 
The variability period, or a provisional estimate of such a period, derived on the basis of the Hipparcos data (possibly in combination with ground-based observations) and expressed in days, with a precision of 0.01 days.

HIP_Var_Type: 
The variability type: the sources of scatter in the photometric data are various, and this flag indicates the origin of the extra scatter, which may be astrophysical, or, in some cases, instrumental. See Section 1.3, Appendix 2 of the published Hipparcos Catalog for a more detailed description. Amongst astrophysical sources of variability, this parameter only distinguishes between 'M' (micro-variables), 'P' (periodic variables), and 'U' (unsolved variables). Further variability details for the periodic or unsolved variables are included in the Variability Annex. The flag takes the following values:

       C : no variability detected ("constant")
       D : duplicity-induced variability
       M : possibly micro-variable, with amplitude < 0.03 mag (stars
           classified with high confidence as micro-variable are flagged U)
       P : periodic variable
       R : the V-I colour index was revised during the variability analysis
       U : unsolved variable which does not fall in the other categories;
           this class also includes irregular or semi-regular variables,
           and possibly varaibles with amplitude > or ~ 0.03 mag
         : a blank indicates that the entry could not be classified as
           variable or constant with any degree of certainty
  
Var_Data_Annex: 
A Variability Annex flag indicating the existence of additional tabular data in the Variability Annex, where '1' means that additional data are provided in a table of periodic variables, and '2' means that additional data are provided in a table of 'unsolved' variables.

Var_Curv_Annex: 
A Variability Annex flag indicating the existence of a light curve, or a folded light curve, in the Variability Annex, where 'A' means the light curve is folded, and 'B' or 'C' mean that the light curve is NOT folded.

CCDM_ID: 
The Catalog of Components of Double and Multiple Stars (CCDM) identifier.

CCDM_History: 
The historical status of the CCDM identifier. The flag takes the following values:

       H : system determined as double or multiple by the Hipparcos
           observations, and was previously unknown as double or multiple
       I : system previously identified as multiple, as given in Annex 1
           of the Hipparcos Input Catalog (HIC)
       M : miscellaneous (system had been previously identified, after
           publication of the HIC, using other more recently available
           catalogs and compilations)
  
CCDM_N_Entries: 
The number of separate catalog entries with the same CCDM identifier.

CCDM_N_Comp: 
The number of components into which the entry was resolved as a result of the satellite observations and data reductions.

Dbl_Mult_Annex: 
The Double and Multiple Systems Annex flag. This indicates that further details of this system are given in one of the 5 (mutually exclusive) parts of the Double and Multiple Systems Annex labelled as follows:

       C : solutions for the components
       G : acceleration or higher order terms
       O : orbital solutions
       V : variability-induced movers (apparent motion arises from variability
           of one of the components of a double system)
       X : stochastic solution (probably astrometric binaries of short period)
  
Astrom_Mult_Source: 
A flag for the source of the absolute astrometry. This parameter qualifies the source of the astrometric parameters for some of the entries with a value of 'C' for the parameter Dbl_Mult_Annex. The values are as follows:

       P : primary target of a 2- or 3-pointing system
       F : secondary or tertiary of a 2- or 3-pointing 'fixed' system
           (common parallax and proper motions)
       I : secondary or tertiary of a 2- or 3-pointing 'independent'
           system (no constraints on parallax or proper motions)
       L : secondary or tertiary of a 2- or 3-pointing 'linear' system
           (common parallax)
       S : astrometric parameters from 'single-star merging' process.
  
Dbl_Soln_Qual: 
A solution quality flag which indicates the reliability of the double or multiple star solution, and is set for all entries in Part C of the Double and Multiple Systems Annex. The flags can be understood as follows:

       A: 'good', or reliable solution
       B: 'fair', or moderately reliable solution
       C: 'poor', or less reliable solution
       D: uncertain solution
       S: suspected non-single, i.e., possible double or multiple,
          although no significant or convincing non-single star solution
          was found
  
Dbl_Ref_ID: 
Component designation for the double star parameters, Dbl_theta, dbl_rho, etc. The first letter gives the 'reference' component, and the second letter gives the subsidiary component. In the case of the Hipparcos observations, the reference component is always defined to be the brighter component (in median H_P) such that the magnitude difference between the components (Diff_Hip_Mag) is always positive.

Dbl_Theta: 
The rounded value for the position angle between the components specified in the Dbl_Ref_id field, expressed in degrees (in the usual sense measured counterclockwise from North).

Dbl_Rho: 
The rounded value for the angular separation between the components specified in the Dbl_Ref_id field, expressed in arcseconds.

Rho_Error: 
The standard error of the angular separation, Dbl_Rho, given in arcseconds.

Diff_Hip_Mag: 
The Hipparcos magnitude difference of the components specified in the Dbl_Ref_id field, expressed in magnitudes.

DHip_Mag_Error: 
The standard error of the Hipparcos magnitude difference, expressed in magnitudes.

Survey_Star: 
A flag indicating a `survey' star. The `survey' was the basic list of bright stars added to and merged with the total list of proposed stars, to provide a stellar sample (almost) complete to well-defined limits. A flag 'S' indicates that the entry is contained within this `survey', whose limiting magnitude is a function of the stars's spectral type and galactic latitude b and is defined by:

     V <= 7.9 + 1.1 x |sin b| for spectral types earlier or equal to G5
     V <= 7.3 + 1.1 x |sin b| for spectral types later than G5
  
If no spectral data were available, the break was taken at (B-V) = 0.8 mag.
ID_Chart
A flag indicating an identification chart. Where identification of a star using ground-based telescopes might prove difficult or ambiguous, identification chrats were constructed and are available in Volume 13 of the printed catalog. Charts correspond to the object observed by the satellite (i.e., at the posotion given in this catalog), even if it was not the intended target. The flag takes the following values: 'D' for charts produced directly from the STScI Digitized Sky Survey (776 entries) or 'G' for charts constructed from the Guide Star Catalog (10877 entries).

Notes
A flag indicating a note is given at the end of the volume(s) in the printed catalog. The flag has the following meaning:

       D : double and multiple systems note only (Volume 10)
       G : general note only (Volumes 5-9)
       P : photometric (including variability) notes only (Volume 11)
       W : D + P only
       X : D + G only
       Y : G + P only
       Z : D + G + P
  
HD_ID: 
HD/HDE/HDEC identifier (CDS Catalog <III 135>).

BD_ID: 
Bonner Durchmusterung (BD) identifier (CDS Catalogs <I 119>, <I 122>). BD identifiers, unlike the CoD and CPD identifiers, may carry a suffix letter for additional stars, i.e., stars with suffixes 'A', "B', 'P', or 'S': these stars were added to the BD Catalog after the original numbering was made, and such suffixes do not imply that the entry is a component of a double or multiple system.

CoD_ID: 
Cordoba Durchmusterung (CoD) identifier (CDS Catalog <I 114>).

CPD_ID: 
Cape Photographic Durchmusterung (CPD) identifier (CDS Catalog <I 108>).

VI_Color_Reduct: 
The (V-I) color index used for the photometric processing (not necessarily the same as the `final' value given in the parameter VI_mag).

Spect_Type: 
The MK or HD spectral type acquired from ground-based compilations and primarily taken from the Hipparcos Input Catalog, with some updates, especially for variable stars.

Spect_Type_Source: 
The source of the spectral type. The flag indicates the source as follows:

   1 : Michigan catalogue for the HD stars, vol. 1 (Houk+, 1975) <III/31>
   2 : Michigan catalogue for the HD stars, vol. 2 (Houk, 1978) <III/51>
   3 : Michigan Catalogue for the HD stars, vol. 3 (Houk, 1982) <III/80>
   4 : Michigan Catalogue for the HD stars, vol. 4 (Houk+, 1988) <III/133>
   G : updated after publication of the HIC <I/196>
   K : General Catalog of Variable Stars, 4th Ed. (Kholopov+ 1988) <II/139>
   S : SIMBAD database at http://cdsweb.u-strasbg.fr/Simbad.html
   X : Miscellaneous
     : A blank entry has no corresponding information.
  
Class: 
The Browse classification created by the HEASARC based on the value of the spect_type parameter.
