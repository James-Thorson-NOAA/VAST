
#' @export
Print_Message = function( whichmessage ){
  if( whichmessage=="GSL_american_plaice" ){
  message(
    "This data set contains data for American plaice in the Gulf of St. Lawrence, as provided by Hugues Benoit.
    The data are from the Gulf Region September Bottom-Trawl survey Database, from the Aquatic Resources Division, Science Branch, Department of Fisheries and Oceans
    "
  )
  }
  if( whichmessage=="GB_haddock" ){
  message(
    "This haddock data set is from the Northeast Fishery Science Center's Spring and Fall Bottom Trawl Survey.

    The Georges Bank haddock data comprise strata 13-25, 29, 30.
    Data were extracted based on SHG criteria (136) for years <=2008, and for TOGA criteria (132x) for years >=2009.

    Door and vessel calibrations were applied where appropriate during the Albatross survey years.  These calibrations are:
    Doors: 1.510 for catch in weight; 1.490 for catch in numbers;
    Vessel: 0.790 for catch in weight; 0.820 for catch in numbers;
    Gear: no conversion was applied (none available).

    In years>=2009, length-season based calibration factors were applied to account for the transition from the Albatross IV to the Henry Bigelow Research Survey Vessel. These calibrations are applied to numbers at length:
    length	calibration_factor
    1	2.626168725
    2	2.626168725
    3	2.626168725
    4	2.626168725
    5	2.626168725
    6	2.626168725
    7	2.626168725
    8	2.626168725
    9	2.626168725
    10	2.626168725
    11	2.626168725
    12	2.626168725
    13	2.626168725
    14	2.626168725
    15	2.626168725
    16	2.626168725
    17	2.626168725
    18	2.626168725
    19	2.580550826
    20	2.534932926
    21	2.489315027
    22	2.443697127
    23	2.398079228
    24	2.352461328
    25	2.306843429
    26	2.261225529
    27	2.215607629
    28	2.16998973
    29	2.12437183
    30	2.078753931
    31	2.033136031
    32	1.987518132
    33	1.941900232
    34	1.896282333
    35	1.850664433
    36	1.805046534
    37	1.759428634
    38	1.713810734
    39	1.668192835
    40	1.622574935
    41	1.576957036
    42	1.531339136
    43	1.485721237
    44	1.440103337
    45	1.394485438
    46	1.348867538
    47	1.303249639
    48	1.257631739
    49	1.21201384
    50	1.16639594
    >=51	1.163989794



    Records where catch in number for haddock were >0 but had corresponding weight =0 were dropped. Also, if catch weight was >0 but length was 0, that record was dropped.  Such occurrences are noted by season.
    *Fall
    240 length records (148 distinct tows) had catchwt=0;
    the length range is 6-42 cm (most are 22 cm or below)

    *Spring
    one station had length=0 but catchwt>0:
    Cruise	Strata	Tow	Station	SPP	SEX	EXPFAC	CATCHWT	CATCHNUM	LENGTH	NUMLEN
    199902	1160	1	258	74	0	1	0.9	1	0	1

    29 length records (20 distinct tows) had catchwt=0;
    the length range 3-86 cm (most are 25 cm or below)


    NOAA Northeast Fisheries Science Center (NEFSC) Database Attributes for raw Haddock Data from Fall 1963 (Spring 1968) to Fall 2014 (Spring 2015)

    COMMON_NAME
    Species common name as assigned in the NEFSC fisheries independent data system.

    SVSPP
    A three digit species code as assigned in the NEFSC fisheries independent database system.

    PURPOSE_CODE
    A two digit code identifying the overall purpose of survey cruise types as they pertain to qualifying data storage and retrieval

    SEASON
    A six character field defining the season a survey cruise type was conducted.  NOTE:  Not all performed survey work is conducted
    	within the dates of standard seasonal definition.  Overlap can occur, but most work is performed within the season defined for a survey.

    CRUISE
    A 6 digit alphanumeric field identifying a unique cruise.

    STRATA
    A five digit alphanumeric code identifying the stratified random sampling areas used in the NEFSC standard fisheries independent trawl surveys.

    TOW
    A three digit alphanumeric field identifying a completed survey trawling operation within a specified stratum.

    STATION
    A four digit alphanumeric field representing the sequential order in which randomly chosen survey trawling operations have been completed.

    SHG
    A three digit alphanumeric field derived from the concatenation of three separate fields used to provide a subjective representation of the
    	quality of a tow.  The S, or station type is used to define the type of tow being performed (e.g., stratified random, non-random, comparison haul, etc.).
    	The H, or haul value represents the relative success of a survey operation (e.g., representative tow, non-representative, may or may not be representative
    	due to gear problems or tow duration). The G, or gear condition identifies possible problems due to gear damage or obstructions in the gear. A combined value
    	equal to 136 or less is the standard performance value used to include survey trawling operations in an analysis. Tows greater than 136 are considered non-representative
    	and are not used in standard analysis.

    TOGA
    A four digit alphanumeric code used to replace the SHG performance indicator upon the replacement of the RV Albatross IV with the FSV Henry B. Bigelow.
    	TOGA is a detailed analysis of survey trawl and vessel performance during each tow, utilizing available data from trawl mensuration systems and vessel sensors
    	routinely logged by the Scientific Computing System (SCS). Tolerance limits and optimal values were calculated from data collected during the NEFSC calibration
    	experiments. These tolerance limits are intended to promote consistency of trawl geometry and towing procedure to validate comparison of the collected trawl survey
    	data with results from the calibration experiments.

    TOWDATE
    A date field representing the UTC date and time the current survey trawling operation has started sampling.

    YEAR
    A four digit alphanumeric field extracted from the TOWDATE field representing the UTC year when a survey trawling operation was conducted.

    MON
    A two digit alphanumeric field extracted from the TOWDATE field representing the UTC month when a survey trawling operation was conducted.

    DAY
    A two digit alphanumeric field extracted from the TOWDATE field representing the UTC day when a survey trawling operation was conducted.

    TIME
    An eight alphanumeric field extracted from the TOWDATE field representing the UTC time when a survey trawling operation was conducted.

    DEPTH
    A four digit number recording the average depth, to the nearest meter, during a survey trawling operation. For surveys starting with the spring of 2001,
    	the depth is recorded every 10 seconds and averaged upon completion of a survey trawling operation. For surveys prior to spring 2001,   avgdepth is derived
    	using the fields SETDEPTH and ENDDEPTH (E.g., (SETDEPTH+ENDDEPTH)/2 avgdepth). When the setdepth is null, then the enddepth is used as the avgdepth and vice
    	versa. If neither setdepth or enddepth is recorded, then the average depth is recorded as null.

    BOT_TEMP
    A three digit numeric field measuring water temperature to the nearest tenth of a degree Celsius near the bottom.  Temperature was recorded using a mechanical
    	bathythermograph in the early years; an expendable bathythermograph (XBT) in the middle years; a CTD from the early 1990s to present.

    CATCHWT
    A nine digit field measuring the aggregate weight of a species.  Catch weights were recorded using a beam balance until the early 1990s measured to the nearest
    	tenth of a kg.  A Marel motion compensating scale has been used to capture weights starting in the early 1990s to present.  The weights included in data files
    	have been adjusted to account for gear and vessel changes over the entire time series where noted.

    CATCHNUM
    A nine digit field identifying the total number of actual or expanded (sub-sampling) fish captured during a trawl operation.  The numbers included in this file
    	have been adjusted to account for gear and vessel changes over the entire time series where noted.

    LATITUDE
    A ten digit number representing degrees decimal minutes of latitude.  Precision of this value has changed over time based on changes in technology.

    LONGITUDE
    A ten digit number representing degrees decimal minutes of longitude.  Precision of this value has changed over time based on changes in technology
    "
    )

  }
}
