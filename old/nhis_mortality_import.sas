*********************************************************************
FEBRUARY 16, 2010  9:34 AM
 
 THIS IS AN EXAMPLE OF A PROGRAM THAT CAN BE USED TO READ IN THE 
 NHIS 1997 PUBLIC-USE LINKED MORTALITY ASCII FILE AND RUN 
 FREQUENCIES ON THE DATA.
 
 NOTE: THE FORMAT DEFINITIONS GIVEN BELOW WILL RESULT IN
       PROCEDURE OUTPUT SHOWING VALUES THAT HAVE BEEN
       GROUPED AS THEY ARE SHOWN IN THE FILE LAYOUT
       DOCUMENTATION
 
NHIS Linked Mortality--2010

Total Records, NHIS 1997: 103,477
***********************************************************************

TO DOWNLOAD AND SAVE THE NHIS 1997 PUBLIC-USE LINKED MORTALITY FILE TO
YOUR HARD DRIVE, FOLLOW THESE STEPS:

*STEP 1: DESIGNATE A FOLDER ON YOUR HARD DRIVE TO DOWNLOAD THE NHIS 1997
         PUBLIC USE LINKED MORTALITY FILE. IN THIS EXAMPLE, THE DATA WILL
         BE SAVED TO: 'C:\Public Use Data\NHIS\DATA'

*STEP 2: GO TO THE NCHS WEB SITE (LINK BELOW) AND DOWNLOAD THE NHIS 1997
         PUBLIC-USE LINKED MORTALITY ASCII FILE TO THE FOLDER 
         'C:\Public Use Data\NHIS\DATA'

ftp://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/NHIS97_MORT_PUBLIC_USE_2010.DAT

***********************************************************************;

libname sk "/home/kovalchs/nhis/mortality/";

* DEFINE VARIABLE VALUES FOR REPORTS;

PROC FORMAT;

  VALUE ELIGFMT
    1 = "Eligible"
    2 = "Under age 18"
    3 = "Ineligible" ;

  VALUE MORTFMT
    0 = "Assumed alive"
    1 = "Assumed deceased"
    . = "Ineligible or under age 18";

  VALUE MRSRCFMT
  	0 = "No"
	1 = "Yes"
	9 = "Not Applicable"
	. = "Ineligible, under age 18 or assumed alive";

 VALUE CAUSEFMT
  	0 = "No"
	1 = "Yes"
	. = "Ineligible, under age 18 or assumed alive";

  VALUE QRTFMT
    1 = "January - March"
    2 = "April   - June"
    3 = "July    - September"
    4 = "October - December" 
    . = "Ineligible, under age 18 or assumed alive";

  VALUE DODYFMT
    . = "Ineligible, under age 18 or assumed alive";

  VALUE FLAGFMT
    0 = "No"
    1 = "Yes"  
    . = "Ineligible, under age 18, assumed alive or no cause data";

RUN ;




DATA sk.NHIS86(drop=HYPERTEN DIABETES HIPFRACT);

	INFILE "/home/kovalchs/nhis/mortality/NHIS86_MORT_PUBLIC_USE_2010.DAT"  LRECL = 61 PAD MISSOVER ;

   * INPUT ALL VARIABLES;
   INPUT
     	PUBLICID 		$  1 - 14
    	PUBLICID_2 		$ 15 - 28   /* THIS IS A BLANK VARIABLE FOR NHIS 1986-1996, 2004 */
    	ELIGSTAT    	29
    	MORTSTAT   		30
		MORTSRCE_NDI	31
		MORTSRCE_SSA	32
		MORTSRCE_CMS	33
     	DODQTR      	34
     	DODYEAR     	35 - 38
		CAUSEAVL		39
     	UCOD_113		$ 40 - 42
     	DIABETES    	43
     	HYPERTEN    	44
     	HIPFRACT    	45
		WGT_NEW    		46 - 53
		SA_WGT_NEW		54 - 61    /* THIS IS A BLANK VARIABLE FOR NHIS 1986-1996 */
     ;



   * DEFINE VARIABLE LABELS;
   LABEL
		PUBLICID		= 'Public-use ID'
		ELIGSTAT		= 'Eligibility Status for Mortality Follow-up'
		MORTSTAT 		= 'Final Mortality Status' 
		MORTSRCE_NDI	= "Mortality Source: NDI Match"
		MORTSRCE_SSA	= "Mortality Source: SSA Information"
		MORTSRCE_CMS	= "Mortality Source: CMS Information"
		DODQTR			= "Quarter of Death"
		DODYEAR			= "Year of Death"
		CAUSEAVL		= "Cause of Death Data Available"
		UCOD_113		= 'Underlying Cause of Death 113 Groups All Years (ICD-10)'
 		DIABETES		= 'Diabetes Flag from Multiple Cause of Death (MCOD)'
		HYPERTEN		= 'Hypertension Flag from Multiple Cause of Death (MCOD)'
		HIPFRACT		= 'Hip Fracture Flag from Multiple Cause of Death (MCOD)'
		WGT_NEW			= "Weight Adjusted for Ineligible Respondents Person-level Sample Weight"
		SA_WGT_NEW		= "Weight Adjusted for Ineligible Respondents Sample Adult Sample Weight"

     ;

   * ASSOCIATE VARIABLES WITH FORMAT VALUES;
   FORMAT    
		ELIGSTAT 		ELIGFMT.          
		MORTSTAT 		MORTFMT.
     	MORTSRCE_NDI 	MRSRCFMT.
     	MORTSRCE_SSA 	MRSRCFMT.
     	MORTSRCE_CMS 	MRSRCFMT.
		CAUSEAVL 		CAUSEFMT.
     	DODQTR   		QRTFMT.           
		DODYEAR  		DODYFMT.
     	DIABETES 		FLAGFMT.          
		HYPERTEN 		FLAGFMT. 
     	HIPFRACT 		FLAGFMT. ;
RUN;

DATA sk.NHIS87(drop=HYPERTEN DIABETES HIPFRACT);

	INFILE "/home/kovalchs/nhis/mortality/NHIS87_MORT_PUBLIC_USE_2010.DAT"  LRECL = 61 PAD MISSOVER ;

   * INPUT ALL VARIABLES;
   INPUT
     	PUBLICID 		$  1 - 14
    	PUBLICID_2 		$ 15 - 28   /* THIS IS A BLANK VARIABLE FOR NHIS 1986-1996, 2004 */
    	ELIGSTAT    	29
    	MORTSTAT   		30
		MORTSRCE_NDI	31
		MORTSRCE_SSA	32
		MORTSRCE_CMS	33
     	DODQTR      	34
     	DODYEAR     	35 - 38
		CAUSEAVL		39
     	UCOD_113		$ 40 - 42
     	DIABETES    	43
     	HYPERTEN    	44
     	HIPFRACT    	45
		WGT_NEW    		46 - 53
		SA_WGT_NEW		54 - 61    /* THIS IS A BLANK VARIABLE FOR NHIS 1986-1996 */
     ;



   * DEFINE VARIABLE LABELS;
   LABEL
		PUBLICID		= 'Public-use ID'
		ELIGSTAT		= 'Eligibility Status for Mortality Follow-up'
		MORTSTAT 		= 'Final Mortality Status' 
		MORTSRCE_NDI	= "Mortality Source: NDI Match"
		MORTSRCE_SSA	= "Mortality Source: SSA Information"
		MORTSRCE_CMS	= "Mortality Source: CMS Information"
		DODQTR			= "Quarter of Death"
		DODYEAR			= "Year of Death"
		CAUSEAVL		= "Cause of Death Data Available"
		UCOD_113		= 'Underlying Cause of Death 113 Groups All Years (ICD-10)'
 		DIABETES		= 'Diabetes Flag from Multiple Cause of Death (MCOD)'
		HYPERTEN		= 'Hypertension Flag from Multiple Cause of Death (MCOD)'
		HIPFRACT		= 'Hip Fracture Flag from Multiple Cause of Death (MCOD)'
		WGT_NEW			= "Weight Adjusted for Ineligible Respondents Person-level Sample Weight"
		SA_WGT_NEW		= "Weight Adjusted for Ineligible Respondents Sample Adult Sample Weight"

     ;

   * ASSOCIATE VARIABLES WITH FORMAT VALUES;
   FORMAT    
		ELIGSTAT 		ELIGFMT.          
		MORTSTAT 		MORTFMT.
     	MORTSRCE_NDI 	MRSRCFMT.
     	MORTSRCE_SSA 	MRSRCFMT.
     	MORTSRCE_CMS 	MRSRCFMT.
		CAUSEAVL 		CAUSEFMT.
     	DODQTR   		QRTFMT.           
		DODYEAR  		DODYFMT.
     	DIABETES 		FLAGFMT.          
		HYPERTEN 		FLAGFMT. 
     	HIPFRACT 		FLAGFMT. ;
RUN;

DATA sk.NHIS88(drop=HYPERTEN DIABETES HIPFRACT);

	INFILE "/home/kovalchs/nhis/mortality/NHIS88_MORT_PUBLIC_USE_2010.DAT"  LRECL = 61 PAD MISSOVER ;

   * INPUT ALL VARIABLES;
   INPUT
     	PUBLICID 		$  1 - 14
    	PUBLICID_2 		$ 15 - 28   /* THIS IS A BLANK VARIABLE FOR NHIS 1986-1996, 2004 */
    	ELIGSTAT    	29
    	MORTSTAT   		30
		MORTSRCE_NDI	31
		MORTSRCE_SSA	32
		MORTSRCE_CMS	33
     	DODQTR      	34
     	DODYEAR     	35 - 38
		CAUSEAVL		39
     	UCOD_113		$ 40 - 42
     	DIABETES    	43
     	HYPERTEN    	44
     	HIPFRACT    	45
		WGT_NEW    		46 - 53
		SA_WGT_NEW		54 - 61    /* THIS IS A BLANK VARIABLE FOR NHIS 1986-1996 */
     ;



   * DEFINE VARIABLE LABELS;
   LABEL
		PUBLICID		= 'Public-use ID'
		ELIGSTAT		= 'Eligibility Status for Mortality Follow-up'
		MORTSTAT 		= 'Final Mortality Status' 
		MORTSRCE_NDI	= "Mortality Source: NDI Match"
		MORTSRCE_SSA	= "Mortality Source: SSA Information"
		MORTSRCE_CMS	= "Mortality Source: CMS Information"
		DODQTR			= "Quarter of Death"
		DODYEAR			= "Year of Death"
		CAUSEAVL		= "Cause of Death Data Available"
		UCOD_113		= 'Underlying Cause of Death 113 Groups All Years (ICD-10)'
 		DIABETES		= 'Diabetes Flag from Multiple Cause of Death (MCOD)'
		HYPERTEN		= 'Hypertension Flag from Multiple Cause of Death (MCOD)'
		HIPFRACT		= 'Hip Fracture Flag from Multiple Cause of Death (MCOD)'
		WGT_NEW			= "Weight Adjusted for Ineligible Respondents Person-level Sample Weight"
		SA_WGT_NEW		= "Weight Adjusted for Ineligible Respondents Sample Adult Sample Weight"

     ;

   * ASSOCIATE VARIABLES WITH FORMAT VALUES;
   FORMAT    
		ELIGSTAT 		ELIGFMT.          
		MORTSTAT 		MORTFMT.
     	MORTSRCE_NDI 	MRSRCFMT.
     	MORTSRCE_SSA 	MRSRCFMT.
     	MORTSRCE_CMS 	MRSRCFMT.
		CAUSEAVL 		CAUSEFMT.
     	DODQTR   		QRTFMT.           
		DODYEAR  		DODYFMT.
     	DIABETES 		FLAGFMT.          
		HYPERTEN 		FLAGFMT. 
     	HIPFRACT 		FLAGFMT. ;
RUN;

DATA sk.NHIS89(drop=HYPERTEN DIABETES HIPFRACT);

	INFILE "/home/kovalchs/nhis/mortality/NHIS89_MORT_PUBLIC_USE_2010.DAT"  LRECL = 61 PAD MISSOVER ;

   * INPUT ALL VARIABLES;
   INPUT
     	PUBLICID 		$  1 - 14
    	PUBLICID_2 		$ 15 - 28   /* THIS IS A BLANK VARIABLE FOR NHIS 1986-1996, 2004 */
    	ELIGSTAT    	29
    	MORTSTAT   		30
		MORTSRCE_NDI	31
		MORTSRCE_SSA	32
		MORTSRCE_CMS	33
     	DODQTR      	34
     	DODYEAR     	35 - 38
		CAUSEAVL		39
     	UCOD_113		$ 40 - 42
     	DIABETES    	43
     	HYPERTEN    	44
     	HIPFRACT    	45
		WGT_NEW    		46 - 53
		SA_WGT_NEW		54 - 61    /* THIS IS A BLANK VARIABLE FOR NHIS 1986-1996 */
     ;



   * DEFINE VARIABLE LABELS;
   LABEL
		PUBLICID		= 'Public-use ID'
		ELIGSTAT		= 'Eligibility Status for Mortality Follow-up'
		MORTSTAT 		= 'Final Mortality Status' 
		MORTSRCE_NDI	= "Mortality Source: NDI Match"
		MORTSRCE_SSA	= "Mortality Source: SSA Information"
		MORTSRCE_CMS	= "Mortality Source: CMS Information"
		DODQTR			= "Quarter of Death"
		DODYEAR			= "Year of Death"
		CAUSEAVL		= "Cause of Death Data Available"
		UCOD_113		= 'Underlying Cause of Death 113 Groups All Years (ICD-10)'
 		DIABETES		= 'Diabetes Flag from Multiple Cause of Death (MCOD)'
		HYPERTEN		= 'Hypertension Flag from Multiple Cause of Death (MCOD)'
		HIPFRACT		= 'Hip Fracture Flag from Multiple Cause of Death (MCOD)'
		WGT_NEW			= "Weight Adjusted for Ineligible Respondents Person-level Sample Weight"
		SA_WGT_NEW		= "Weight Adjusted for Ineligible Respondents Sample Adult Sample Weight"

     ;

   * ASSOCIATE VARIABLES WITH FORMAT VALUES;
   FORMAT    
		ELIGSTAT 		ELIGFMT.          
		MORTSTAT 		MORTFMT.
     	MORTSRCE_NDI 	MRSRCFMT.
     	MORTSRCE_SSA 	MRSRCFMT.
     	MORTSRCE_CMS 	MRSRCFMT.
		CAUSEAVL 		CAUSEFMT.
     	DODQTR   		QRTFMT.           
		DODYEAR  		DODYFMT.
     	DIABETES 		FLAGFMT.          
		HYPERTEN 		FLAGFMT. 
     	HIPFRACT 		FLAGFMT. ;
RUN;

DATA sk.NHIS90(drop=HYPERTEN DIABETES HIPFRACT);

	INFILE "/home/kovalchs/nhis/mortality/NHIS90_MORT_PUBLIC_USE_2010.DAT"  LRECL = 61 PAD MISSOVER ;

   * INPUT ALL VARIABLES;
   INPUT
     	PUBLICID 		$  1 - 14
    	PUBLICID_2 		$ 15 - 28   /* THIS IS A BLANK VARIABLE FOR NHIS 1986-1996, 2004 */
    	ELIGSTAT    	29
    	MORTSTAT   		30
		MORTSRCE_NDI	31
		MORTSRCE_SSA	32
		MORTSRCE_CMS	33
     	DODQTR      	34
     	DODYEAR     	35 - 38
		CAUSEAVL		39
     	UCOD_113		$ 40 - 42
     	DIABETES    	43
     	HYPERTEN    	44
     	HIPFRACT    	45
		WGT_NEW    		46 - 53
		SA_WGT_NEW		54 - 61    /* THIS IS A BLANK VARIABLE FOR NHIS 1986-1996 */
     ;



   * DEFINE VARIABLE LABELS;
   LABEL
		PUBLICID		= 'Public-use ID'
		ELIGSTAT		= 'Eligibility Status for Mortality Follow-up'
		MORTSTAT 		= 'Final Mortality Status' 
		MORTSRCE_NDI	= "Mortality Source: NDI Match"
		MORTSRCE_SSA	= "Mortality Source: SSA Information"
		MORTSRCE_CMS	= "Mortality Source: CMS Information"
		DODQTR			= "Quarter of Death"
		DODYEAR			= "Year of Death"
		CAUSEAVL		= "Cause of Death Data Available"
		UCOD_113		= 'Underlying Cause of Death 113 Groups All Years (ICD-10)'
 		DIABETES		= 'Diabetes Flag from Multiple Cause of Death (MCOD)'
		HYPERTEN		= 'Hypertension Flag from Multiple Cause of Death (MCOD)'
		HIPFRACT		= 'Hip Fracture Flag from Multiple Cause of Death (MCOD)'
		WGT_NEW			= "Weight Adjusted for Ineligible Respondents Person-level Sample Weight"
		SA_WGT_NEW		= "Weight Adjusted for Ineligible Respondents Sample Adult Sample Weight"

     ;

   * ASSOCIATE VARIABLES WITH FORMAT VALUES;
   FORMAT    
		ELIGSTAT 		ELIGFMT.          
		MORTSTAT 		MORTFMT.
     	MORTSRCE_NDI 	MRSRCFMT.
     	MORTSRCE_SSA 	MRSRCFMT.
     	MORTSRCE_CMS 	MRSRCFMT.
		CAUSEAVL 		CAUSEFMT.
     	DODQTR   		QRTFMT.           
		DODYEAR  		DODYFMT.
     	DIABETES 		FLAGFMT.          
		HYPERTEN 		FLAGFMT. 
     	HIPFRACT 		FLAGFMT. ;
RUN;

DATA sk.NHIS91(drop=HYPERTEN DIABETES HIPFRACT);

	INFILE "/home/kovalchs/nhis/mortality/NHIS91_MORT_PUBLIC_USE_2010.DAT"  LRECL = 61 PAD MISSOVER ;

   * INPUT ALL VARIABLES;
   INPUT
     	PUBLICID 		$  1 - 14
    	PUBLICID_2 		$ 15 - 28   /* THIS IS A BLANK VARIABLE FOR NHIS 1986-1996, 2004 */
    	ELIGSTAT    	29
    	MORTSTAT   		30
		MORTSRCE_NDI	31
		MORTSRCE_SSA	32
		MORTSRCE_CMS	33
     	DODQTR      	34
     	DODYEAR     	35 - 38
		CAUSEAVL		39
     	UCOD_113		$ 40 - 42
     	DIABETES    	43
     	HYPERTEN    	44
     	HIPFRACT    	45
		WGT_NEW    		46 - 53
		SA_WGT_NEW		54 - 61    /* THIS IS A BLANK VARIABLE FOR NHIS 1986-1996 */
     ;



   * DEFINE VARIABLE LABELS;
   LABEL
		PUBLICID		= 'Public-use ID'
		ELIGSTAT		= 'Eligibility Status for Mortality Follow-up'
		MORTSTAT 		= 'Final Mortality Status' 
		MORTSRCE_NDI	= "Mortality Source: NDI Match"
		MORTSRCE_SSA	= "Mortality Source: SSA Information"
		MORTSRCE_CMS	= "Mortality Source: CMS Information"
		DODQTR			= "Quarter of Death"
		DODYEAR			= "Year of Death"
		CAUSEAVL		= "Cause of Death Data Available"
		UCOD_113		= 'Underlying Cause of Death 113 Groups All Years (ICD-10)'
 		DIABETES		= 'Diabetes Flag from Multiple Cause of Death (MCOD)'
		HYPERTEN		= 'Hypertension Flag from Multiple Cause of Death (MCOD)'
		HIPFRACT		= 'Hip Fracture Flag from Multiple Cause of Death (MCOD)'
		WGT_NEW			= "Weight Adjusted for Ineligible Respondents Person-level Sample Weight"
		SA_WGT_NEW		= "Weight Adjusted for Ineligible Respondents Sample Adult Sample Weight"

     ;

   * ASSOCIATE VARIABLES WITH FORMAT VALUES;
   FORMAT    
		ELIGSTAT 		ELIGFMT.          
		MORTSTAT 		MORTFMT.
     	MORTSRCE_NDI 	MRSRCFMT.
     	MORTSRCE_SSA 	MRSRCFMT.
     	MORTSRCE_CMS 	MRSRCFMT.
		CAUSEAVL 		CAUSEFMT.
     	DODQTR   		QRTFMT.           
		DODYEAR  		DODYFMT.
     	DIABETES 		FLAGFMT.          
		HYPERTEN 		FLAGFMT. 
     	HIPFRACT 		FLAGFMT. ;
RUN;

DATA sk.NHIS92(drop=HYPERTEN DIABETES HIPFRACT);

	INFILE "/home/kovalchs/nhis/mortality/NHIS92_MORT_PUBLIC_USE_2010.DAT"  LRECL = 61 PAD MISSOVER ;

   * INPUT ALL VARIABLES;
   INPUT
     	PUBLICID 		$  1 - 14
    	PUBLICID_2 		$ 15 - 28   /* THIS IS A BLANK VARIABLE FOR NHIS 1986-1996, 2004 */
    	ELIGSTAT    	29
    	MORTSTAT   		30
		MORTSRCE_NDI	31
		MORTSRCE_SSA	32
		MORTSRCE_CMS	33
     	DODQTR      	34
     	DODYEAR     	35 - 38
		CAUSEAVL		39
     	UCOD_113		$ 40 - 42
     	DIABETES    	43
     	HYPERTEN    	44
     	HIPFRACT    	45
		WGT_NEW    		46 - 53
		SA_WGT_NEW		54 - 61    /* THIS IS A BLANK VARIABLE FOR NHIS 1986-1996 */
     ;



   * DEFINE VARIABLE LABELS;
   LABEL
		PUBLICID		= 'Public-use ID'
		ELIGSTAT		= 'Eligibility Status for Mortality Follow-up'
		MORTSTAT 		= 'Final Mortality Status' 
		MORTSRCE_NDI	= "Mortality Source: NDI Match"
		MORTSRCE_SSA	= "Mortality Source: SSA Information"
		MORTSRCE_CMS	= "Mortality Source: CMS Information"
		DODQTR			= "Quarter of Death"
		DODYEAR			= "Year of Death"
		CAUSEAVL		= "Cause of Death Data Available"
		UCOD_113		= 'Underlying Cause of Death 113 Groups All Years (ICD-10)'
 		DIABETES		= 'Diabetes Flag from Multiple Cause of Death (MCOD)'
		HYPERTEN		= 'Hypertension Flag from Multiple Cause of Death (MCOD)'
		HIPFRACT		= 'Hip Fracture Flag from Multiple Cause of Death (MCOD)'
		WGT_NEW			= "Weight Adjusted for Ineligible Respondents Person-level Sample Weight"
		SA_WGT_NEW		= "Weight Adjusted for Ineligible Respondents Sample Adult Sample Weight"

     ;

   * ASSOCIATE VARIABLES WITH FORMAT VALUES;
   FORMAT    
		ELIGSTAT 		ELIGFMT.          
		MORTSTAT 		MORTFMT.
     	MORTSRCE_NDI 	MRSRCFMT.
     	MORTSRCE_SSA 	MRSRCFMT.
     	MORTSRCE_CMS 	MRSRCFMT.
		CAUSEAVL 		CAUSEFMT.
     	DODQTR   		QRTFMT.           
		DODYEAR  		DODYFMT.
     	DIABETES 		FLAGFMT.          
		HYPERTEN 		FLAGFMT. 
     	HIPFRACT 		FLAGFMT. ;
RUN;

DATA sk.NHIS93(drop=HYPERTEN DIABETES HIPFRACT);

	INFILE "/home/kovalchs/nhis/mortality/NHIS93_MORT_PUBLIC_USE_2010.DAT"  LRECL = 61 PAD MISSOVER ;

   * INPUT ALL VARIABLES;
   INPUT
     	PUBLICID 		$  1 - 14
    	PUBLICID_2 		$ 15 - 28   /* THIS IS A BLANK VARIABLE FOR NHIS 1986-1996, 2004 */
    	ELIGSTAT    	29
    	MORTSTAT   		30
		MORTSRCE_NDI	31
		MORTSRCE_SSA	32
		MORTSRCE_CMS	33
     	DODQTR      	34
     	DODYEAR     	35 - 38
		CAUSEAVL		39
     	UCOD_113		$ 40 - 42
     	DIABETES    	43
     	HYPERTEN    	44
     	HIPFRACT    	45
		WGT_NEW    		46 - 53
		SA_WGT_NEW		54 - 61    /* THIS IS A BLANK VARIABLE FOR NHIS 1986-1996 */
     ;



   * DEFINE VARIABLE LABELS;
   LABEL
		PUBLICID		= 'Public-use ID'
		ELIGSTAT		= 'Eligibility Status for Mortality Follow-up'
		MORTSTAT 		= 'Final Mortality Status' 
		MORTSRCE_NDI	= "Mortality Source: NDI Match"
		MORTSRCE_SSA	= "Mortality Source: SSA Information"
		MORTSRCE_CMS	= "Mortality Source: CMS Information"
		DODQTR			= "Quarter of Death"
		DODYEAR			= "Year of Death"
		CAUSEAVL		= "Cause of Death Data Available"
		UCOD_113		= 'Underlying Cause of Death 113 Groups All Years (ICD-10)'
 		DIABETES		= 'Diabetes Flag from Multiple Cause of Death (MCOD)'
		HYPERTEN		= 'Hypertension Flag from Multiple Cause of Death (MCOD)'
		HIPFRACT		= 'Hip Fracture Flag from Multiple Cause of Death (MCOD)'
		WGT_NEW			= "Weight Adjusted for Ineligible Respondents Person-level Sample Weight"
		SA_WGT_NEW		= "Weight Adjusted for Ineligible Respondents Sample Adult Sample Weight"

     ;

   * ASSOCIATE VARIABLES WITH FORMAT VALUES;
   FORMAT    
		ELIGSTAT 		ELIGFMT.          
		MORTSTAT 		MORTFMT.
     	MORTSRCE_NDI 	MRSRCFMT.
     	MORTSRCE_SSA 	MRSRCFMT.
     	MORTSRCE_CMS 	MRSRCFMT.
		CAUSEAVL 		CAUSEFMT.
     	DODQTR   		QRTFMT.           
		DODYEAR  		DODYFMT.
     	DIABETES 		FLAGFMT.          
		HYPERTEN 		FLAGFMT. 
     	HIPFRACT 		FLAGFMT. ;
RUN;

DATA sk.NHIS94(drop=HYPERTEN DIABETES HIPFRACT);

	INFILE "/home/kovalchs/nhis/mortality/NHIS94_MORT_PUBLIC_USE_2010.DAT"  LRECL = 61 PAD MISSOVER ;

   * INPUT ALL VARIABLES;
   INPUT
     	PUBLICID 		$  1 - 14
    	PUBLICID_2 		$ 15 - 28   /* THIS IS A BLANK VARIABLE FOR NHIS 1986-1996, 2004 */
    	ELIGSTAT    	29
    	MORTSTAT   		30
		MORTSRCE_NDI	31
		MORTSRCE_SSA	32
		MORTSRCE_CMS	33
     	DODQTR      	34
     	DODYEAR     	35 - 38
		CAUSEAVL		39
     	UCOD_113		$ 40 - 42
     	DIABETES    	43
     	HYPERTEN    	44
     	HIPFRACT    	45
		WGT_NEW    		46 - 53
		SA_WGT_NEW		54 - 61    /* THIS IS A BLANK VARIABLE FOR NHIS 1986-1996 */
     ;



   * DEFINE VARIABLE LABELS;
   LABEL
		PUBLICID		= 'Public-use ID'
		ELIGSTAT		= 'Eligibility Status for Mortality Follow-up'
		MORTSTAT 		= 'Final Mortality Status' 
		MORTSRCE_NDI	= "Mortality Source: NDI Match"
		MORTSRCE_SSA	= "Mortality Source: SSA Information"
		MORTSRCE_CMS	= "Mortality Source: CMS Information"
		DODQTR			= "Quarter of Death"
		DODYEAR			= "Year of Death"
		CAUSEAVL		= "Cause of Death Data Available"
		UCOD_113		= 'Underlying Cause of Death 113 Groups All Years (ICD-10)'
 		DIABETES		= 'Diabetes Flag from Multiple Cause of Death (MCOD)'
		HYPERTEN		= 'Hypertension Flag from Multiple Cause of Death (MCOD)'
		HIPFRACT		= 'Hip Fracture Flag from Multiple Cause of Death (MCOD)'
		WGT_NEW			= "Weight Adjusted for Ineligible Respondents Person-level Sample Weight"
		SA_WGT_NEW		= "Weight Adjusted for Ineligible Respondents Sample Adult Sample Weight"

     ;

   * ASSOCIATE VARIABLES WITH FORMAT VALUES;
   FORMAT    
		ELIGSTAT 		ELIGFMT.          
		MORTSTAT 		MORTFMT.
     	MORTSRCE_NDI 	MRSRCFMT.
     	MORTSRCE_SSA 	MRSRCFMT.
     	MORTSRCE_CMS 	MRSRCFMT.
		CAUSEAVL 		CAUSEFMT.
     	DODQTR   		QRTFMT.           
		DODYEAR  		DODYFMT.
     	DIABETES 		FLAGFMT.          
		HYPERTEN 		FLAGFMT. 
     	HIPFRACT 		FLAGFMT. ;
RUN;

DATA sk.NHIS95(drop=HYPERTEN DIABETES HIPFRACT);

	INFILE "/home/kovalchs/nhis/mortality/NHIS95_MORT_PUBLIC_USE_2010.DAT"  LRECL = 61 PAD MISSOVER ;

   * INPUT ALL VARIABLES;
   INPUT
     	PUBLICID 		$  1 - 14
    	PUBLICID_2 		$ 15 - 28   /* THIS IS A BLANK VARIABLE FOR NHIS 1986-1996, 2004 */
    	ELIGSTAT    	29
    	MORTSTAT   		30
		MORTSRCE_NDI	31
		MORTSRCE_SSA	32
		MORTSRCE_CMS	33
     	DODQTR      	34
     	DODYEAR     	35 - 38
		CAUSEAVL		39
     	UCOD_113		$ 40 - 42
     	DIABETES    	43
     	HYPERTEN    	44
     	HIPFRACT    	45
		WGT_NEW    		46 - 53
		SA_WGT_NEW		54 - 61    /* THIS IS A BLANK VARIABLE FOR NHIS 1986-1996 */
     ;



   * DEFINE VARIABLE LABELS;
   LABEL
		PUBLICID		= 'Public-use ID'
		ELIGSTAT		= 'Eligibility Status for Mortality Follow-up'
		MORTSTAT 		= 'Final Mortality Status' 
		MORTSRCE_NDI	= "Mortality Source: NDI Match"
		MORTSRCE_SSA	= "Mortality Source: SSA Information"
		MORTSRCE_CMS	= "Mortality Source: CMS Information"
		DODQTR			= "Quarter of Death"
		DODYEAR			= "Year of Death"
		CAUSEAVL		= "Cause of Death Data Available"
		UCOD_113		= 'Underlying Cause of Death 113 Groups All Years (ICD-10)'
 		DIABETES		= 'Diabetes Flag from Multiple Cause of Death (MCOD)'
		HYPERTEN		= 'Hypertension Flag from Multiple Cause of Death (MCOD)'
		HIPFRACT		= 'Hip Fracture Flag from Multiple Cause of Death (MCOD)'
		WGT_NEW			= "Weight Adjusted for Ineligible Respondents Person-level Sample Weight"
		SA_WGT_NEW		= "Weight Adjusted for Ineligible Respondents Sample Adult Sample Weight"

     ;

   * ASSOCIATE VARIABLES WITH FORMAT VALUES;
   FORMAT    
		ELIGSTAT 		ELIGFMT.          
		MORTSTAT 		MORTFMT.
     	MORTSRCE_NDI 	MRSRCFMT.
     	MORTSRCE_SSA 	MRSRCFMT.
     	MORTSRCE_CMS 	MRSRCFMT.
		CAUSEAVL 		CAUSEFMT.
     	DODQTR   		QRTFMT.           
		DODYEAR  		DODYFMT.
     	DIABETES 		FLAGFMT.          
		HYPERTEN 		FLAGFMT. 
     	HIPFRACT 		FLAGFMT. ;
RUN;

DATA sk.NHIS96(drop=HYPERTEN DIABETES HIPFRACT);

	INFILE "/home/kovalchs/nhis/mortality/NHIS96_MORT_PUBLIC_USE_2010.DAT"  LRECL = 61 PAD MISSOVER ;

   * INPUT ALL VARIABLES;
   INPUT
     	PUBLICID 		$  1 - 14
    	PUBLICID_2 		$ 15 - 28   /* THIS IS A BLANK VARIABLE FOR NHIS 1986-1996, 2004 */
    	ELIGSTAT    	29
    	MORTSTAT   		30
		MORTSRCE_NDI	31
		MORTSRCE_SSA	32
		MORTSRCE_CMS	33
     	DODQTR      	34
     	DODYEAR     	35 - 38
		CAUSEAVL		39
     	UCOD_113		$ 40 - 42
     	DIABETES    	43
     	HYPERTEN    	44
     	HIPFRACT    	45
		WGT_NEW    		46 - 53
		SA_WGT_NEW		54 - 61    /* THIS IS A BLANK VARIABLE FOR NHIS 1986-1996 */
     ;



   * DEFINE VARIABLE LABELS;
   LABEL
		PUBLICID		= 'Public-use ID'
		ELIGSTAT		= 'Eligibility Status for Mortality Follow-up'
		MORTSTAT 		= 'Final Mortality Status' 
		MORTSRCE_NDI	= "Mortality Source: NDI Match"
		MORTSRCE_SSA	= "Mortality Source: SSA Information"
		MORTSRCE_CMS	= "Mortality Source: CMS Information"
		DODQTR			= "Quarter of Death"
		DODYEAR			= "Year of Death"
		CAUSEAVL		= "Cause of Death Data Available"
		UCOD_113		= 'Underlying Cause of Death 113 Groups All Years (ICD-10)'
 		DIABETES		= 'Diabetes Flag from Multiple Cause of Death (MCOD)'
		HYPERTEN		= 'Hypertension Flag from Multiple Cause of Death (MCOD)'
		HIPFRACT		= 'Hip Fracture Flag from Multiple Cause of Death (MCOD)'
		WGT_NEW			= "Weight Adjusted for Ineligible Respondents Person-level Sample Weight"
		SA_WGT_NEW		= "Weight Adjusted for Ineligible Respondents Sample Adult Sample Weight"

     ;

   * ASSOCIATE VARIABLES WITH FORMAT VALUES;
   FORMAT    
		ELIGSTAT 		ELIGFMT.          
		MORTSTAT 		MORTFMT.
     	MORTSRCE_NDI 	MRSRCFMT.
     	MORTSRCE_SSA 	MRSRCFMT.
     	MORTSRCE_CMS 	MRSRCFMT.
		CAUSEAVL 		CAUSEFMT.
     	DODQTR   		QRTFMT.           
		DODYEAR  		DODYFMT.
     	DIABETES 		FLAGFMT.          
		HYPERTEN 		FLAGFMT. 
     	HIPFRACT 		FLAGFMT. ;
RUN;

DATA sk.NHIS97(drop=HYPERTEN DIABETES HIPFRACT);

	INFILE "/home/kovalchs/nhis/mortality/NHIS97_MORT_PUBLIC_USE_2010.DAT"  LRECL = 61 PAD MISSOVER ;

   * INPUT ALL VARIABLES;
   INPUT
     	PUBLICID 		$  1 - 14
    	PUBLICID_2 		$ 15 - 28   /* THIS IS A BLANK VARIABLE FOR NHIS 1986-1996, 2004 */
    	ELIGSTAT    	29
    	MORTSTAT   		30
		MORTSRCE_NDI	31
		MORTSRCE_SSA	32
		MORTSRCE_CMS	33
     	DODQTR      	34
     	DODYEAR     	35 - 38
		CAUSEAVL		39
     	UCOD_113		$ 40 - 42
     	DIABETES    	43
     	HYPERTEN    	44
     	HIPFRACT    	45
		WGT_NEW    		46 - 53
		SA_WGT_NEW		54 - 61    /* THIS IS A BLANK VARIABLE FOR NHIS 1986-1996 */
     ;



   * DEFINE VARIABLE LABELS;
   LABEL
		PUBLICID		= 'Public-use ID'
		ELIGSTAT		= 'Eligibility Status for Mortality Follow-up'
		MORTSTAT 		= 'Final Mortality Status' 
		MORTSRCE_NDI	= "Mortality Source: NDI Match"
		MORTSRCE_SSA	= "Mortality Source: SSA Information"
		MORTSRCE_CMS	= "Mortality Source: CMS Information"
		DODQTR			= "Quarter of Death"
		DODYEAR			= "Year of Death"
		CAUSEAVL		= "Cause of Death Data Available"
		UCOD_113		= 'Underlying Cause of Death 113 Groups All Years (ICD-10)'
 		DIABETES		= 'Diabetes Flag from Multiple Cause of Death (MCOD)'
		HYPERTEN		= 'Hypertension Flag from Multiple Cause of Death (MCOD)'
		HIPFRACT		= 'Hip Fracture Flag from Multiple Cause of Death (MCOD)'
		WGT_NEW			= "Weight Adjusted for Ineligible Respondents Person-level Sample Weight"
		SA_WGT_NEW		= "Weight Adjusted for Ineligible Respondents Sample Adult Sample Weight"

     ;

   * ASSOCIATE VARIABLES WITH FORMAT VALUES;
   FORMAT    
		ELIGSTAT 		ELIGFMT.          
		MORTSTAT 		MORTFMT.
     	MORTSRCE_NDI 	MRSRCFMT.
     	MORTSRCE_SSA 	MRSRCFMT.
     	MORTSRCE_CMS 	MRSRCFMT.
		CAUSEAVL 		CAUSEFMT.
     	DODQTR   		QRTFMT.           
		DODYEAR  		DODYFMT.
     	DIABETES 		FLAGFMT.          
		HYPERTEN 		FLAGFMT. 
     	HIPFRACT 		FLAGFMT. ;
RUN;

DATA sk.NHIS98(drop=HYPERTEN DIABETES HIPFRACT);

	INFILE "/home/kovalchs/nhis/mortality/NHIS98_MORT_PUBLIC_USE_2010.DAT"  LRECL = 61 PAD MISSOVER ;

   * INPUT ALL VARIABLES;
   INPUT
     	PUBLICID 		$  1 - 14
    	PUBLICID_2 		$ 15 - 28   /* THIS IS A BLANK VARIABLE FOR NHIS 1986-1996, 2004 */
    	ELIGSTAT    	29
    	MORTSTAT   		30
		MORTSRCE_NDI	31
		MORTSRCE_SSA	32
		MORTSRCE_CMS	33
     	DODQTR      	34
     	DODYEAR     	35 - 38
		CAUSEAVL		39
     	UCOD_113		$ 40 - 42
     	DIABETES    	43
     	HYPERTEN    	44
     	HIPFRACT    	45
		WGT_NEW    		46 - 53
		SA_WGT_NEW		54 - 61    /* THIS IS A BLANK VARIABLE FOR NHIS 1986-1996 */
     ;



   * DEFINE VARIABLE LABELS;
   LABEL
		PUBLICID		= 'Public-use ID'
		ELIGSTAT		= 'Eligibility Status for Mortality Follow-up'
		MORTSTAT 		= 'Final Mortality Status' 
		MORTSRCE_NDI	= "Mortality Source: NDI Match"
		MORTSRCE_SSA	= "Mortality Source: SSA Information"
		MORTSRCE_CMS	= "Mortality Source: CMS Information"
		DODQTR			= "Quarter of Death"
		DODYEAR			= "Year of Death"
		CAUSEAVL		= "Cause of Death Data Available"
		UCOD_113		= 'Underlying Cause of Death 113 Groups All Years (ICD-10)'
 		DIABETES		= 'Diabetes Flag from Multiple Cause of Death (MCOD)'
		HYPERTEN		= 'Hypertension Flag from Multiple Cause of Death (MCOD)'
		HIPFRACT		= 'Hip Fracture Flag from Multiple Cause of Death (MCOD)'
		WGT_NEW			= "Weight Adjusted for Ineligible Respondents Person-level Sample Weight"
		SA_WGT_NEW		= "Weight Adjusted for Ineligible Respondents Sample Adult Sample Weight"

     ;

   * ASSOCIATE VARIABLES WITH FORMAT VALUES;
   FORMAT    
		ELIGSTAT 		ELIGFMT.          
		MORTSTAT 		MORTFMT.
     	MORTSRCE_NDI 	MRSRCFMT.
     	MORTSRCE_SSA 	MRSRCFMT.
     	MORTSRCE_CMS 	MRSRCFMT.
		CAUSEAVL 		CAUSEFMT.
     	DODQTR   		QRTFMT.           
		DODYEAR  		DODYFMT.
     	DIABETES 		FLAGFMT.          
		HYPERTEN 		FLAGFMT. 
     	HIPFRACT 		FLAGFMT. ;
RUN;

DATA sk.NHIS99(drop=HYPERTEN DIABETES HIPFRACT);

	INFILE "/home/kovalchs/nhis/mortality/NHIS99_MORT_PUBLIC_USE_2010.DAT"  LRECL = 61 PAD MISSOVER ;

   * INPUT ALL VARIABLES;
   INPUT
     	PUBLICID 		$  1 - 14
    	PUBLICID_2 		$ 15 - 28   /* THIS IS A BLANK VARIABLE FOR NHIS 1986-1996, 2004 */
    	ELIGSTAT    	29
    	MORTSTAT   		30
		MORTSRCE_NDI	31
		MORTSRCE_SSA	32
		MORTSRCE_CMS	33
     	DODQTR      	34
     	DODYEAR     	35 - 38
		CAUSEAVL		39
     	UCOD_113		$ 40 - 42
     	DIABETES    	43
     	HYPERTEN    	44
     	HIPFRACT    	45
		WGT_NEW    		46 - 53
		SA_WGT_NEW		54 - 61    /* THIS IS A BLANK VARIABLE FOR NHIS 1986-1996 */
     ;



   * DEFINE VARIABLE LABELS;
   LABEL
		PUBLICID		= 'Public-use ID'
		ELIGSTAT		= 'Eligibility Status for Mortality Follow-up'
		MORTSTAT 		= 'Final Mortality Status' 
		MORTSRCE_NDI	= "Mortality Source: NDI Match"
		MORTSRCE_SSA	= "Mortality Source: SSA Information"
		MORTSRCE_CMS	= "Mortality Source: CMS Information"
		DODQTR			= "Quarter of Death"
		DODYEAR			= "Year of Death"
		CAUSEAVL		= "Cause of Death Data Available"
		UCOD_113		= 'Underlying Cause of Death 113 Groups All Years (ICD-10)'
 		DIABETES		= 'Diabetes Flag from Multiple Cause of Death (MCOD)'
		HYPERTEN		= 'Hypertension Flag from Multiple Cause of Death (MCOD)'
		HIPFRACT		= 'Hip Fracture Flag from Multiple Cause of Death (MCOD)'
		WGT_NEW			= "Weight Adjusted for Ineligible Respondents Person-level Sample Weight"
		SA_WGT_NEW		= "Weight Adjusted for Ineligible Respondents Sample Adult Sample Weight"

     ;

   * ASSOCIATE VARIABLES WITH FORMAT VALUES;
   FORMAT    
		ELIGSTAT 		ELIGFMT.          
		MORTSTAT 		MORTFMT.
     	MORTSRCE_NDI 	MRSRCFMT.
     	MORTSRCE_SSA 	MRSRCFMT.
     	MORTSRCE_CMS 	MRSRCFMT.
		CAUSEAVL 		CAUSEFMT.
     	DODQTR   		QRTFMT.           
		DODYEAR  		DODYFMT.
     	DIABETES 		FLAGFMT.          
		HYPERTEN 		FLAGFMT. 
     	HIPFRACT 		FLAGFMT. ;
RUN;

DATA sk.NHIS00(drop=HYPERTEN DIABETES HIPFRACT);

	INFILE "/home/kovalchs/nhis/mortality/NHIS00_MORT_PUBLIC_USE_2010.DAT"  LRECL = 61 PAD MISSOVER ;

   * INPUT ALL VARIABLES;
   INPUT
     	PUBLICID 		$  1 - 14
    	PUBLICID_2 		$ 15 - 28   /* THIS IS A BLANK VARIABLE FOR NHIS 1986-1996, 2004 */
    	ELIGSTAT    	29
    	MORTSTAT   		30
		MORTSRCE_NDI	31
		MORTSRCE_SSA	32
		MORTSRCE_CMS	33
     	DODQTR      	34
     	DODYEAR     	35 - 38
		CAUSEAVL		39
     	UCOD_113		$ 40 - 42
     	DIABETES    	43
     	HYPERTEN    	44
     	HIPFRACT    	45
		WGT_NEW    		46 - 53
		SA_WGT_NEW		54 - 61    /* THIS IS A BLANK VARIABLE FOR NHIS 1986-1996 */
     ;



   * DEFINE VARIABLE LABELS;
   LABEL
		PUBLICID		= 'Public-use ID'
		ELIGSTAT		= 'Eligibility Status for Mortality Follow-up'
		MORTSTAT 		= 'Final Mortality Status' 
		MORTSRCE_NDI	= "Mortality Source: NDI Match"
		MORTSRCE_SSA	= "Mortality Source: SSA Information"
		MORTSRCE_CMS	= "Mortality Source: CMS Information"
		DODQTR			= "Quarter of Death"
		DODYEAR			= "Year of Death"
		CAUSEAVL		= "Cause of Death Data Available"
		UCOD_113		= 'Underlying Cause of Death 113 Groups All Years (ICD-10)'
 		DIABETES		= 'Diabetes Flag from Multiple Cause of Death (MCOD)'
		HYPERTEN		= 'Hypertension Flag from Multiple Cause of Death (MCOD)'
		HIPFRACT		= 'Hip Fracture Flag from Multiple Cause of Death (MCOD)'
		WGT_NEW			= "Weight Adjusted for Ineligible Respondents Person-level Sample Weight"
		SA_WGT_NEW		= "Weight Adjusted for Ineligible Respondents Sample Adult Sample Weight"

     ;

   * ASSOCIATE VARIABLES WITH FORMAT VALUES;
   FORMAT    
		ELIGSTAT 		ELIGFMT.          
		MORTSTAT 		MORTFMT.
     	MORTSRCE_NDI 	MRSRCFMT.
     	MORTSRCE_SSA 	MRSRCFMT.
     	MORTSRCE_CMS 	MRSRCFMT.
		CAUSEAVL 		CAUSEFMT.
     	DODQTR   		QRTFMT.           
		DODYEAR  		DODYFMT.
     	DIABETES 		FLAGFMT.          
		HYPERTEN 		FLAGFMT. 
     	HIPFRACT 		FLAGFMT. ;
RUN;

DATA sk.NHIS01(drop=HYPERTEN DIABETES HIPFRACT);

	INFILE "/home/kovalchs/nhis/mortality/NHIS01_MORT_PUBLIC_USE_2010.DAT"  LRECL = 61 PAD MISSOVER ;

   * INPUT ALL VARIABLES;
   INPUT
     	PUBLICID 		$  1 - 14
    	PUBLICID_2 		$ 15 - 28   /* THIS IS A BLANK VARIABLE FOR NHIS 1986-1996, 2004 */
    	ELIGSTAT    	29
    	MORTSTAT   		30
		MORTSRCE_NDI	31
		MORTSRCE_SSA	32
		MORTSRCE_CMS	33
     	DODQTR      	34
     	DODYEAR     	35 - 38
		CAUSEAVL		39
     	UCOD_113		$ 40 - 42
     	DIABETES    	43
     	HYPERTEN    	44
     	HIPFRACT    	45
		WGT_NEW    		46 - 53
		SA_WGT_NEW		54 - 61    /* THIS IS A BLANK VARIABLE FOR NHIS 1986-1996 */
     ;



   * DEFINE VARIABLE LABELS;
   LABEL
		PUBLICID		= 'Public-use ID'
		ELIGSTAT		= 'Eligibility Status for Mortality Follow-up'
		MORTSTAT 		= 'Final Mortality Status' 
		MORTSRCE_NDI	= "Mortality Source: NDI Match"
		MORTSRCE_SSA	= "Mortality Source: SSA Information"
		MORTSRCE_CMS	= "Mortality Source: CMS Information"
		DODQTR			= "Quarter of Death"
		DODYEAR			= "Year of Death"
		CAUSEAVL		= "Cause of Death Data Available"
		UCOD_113		= 'Underlying Cause of Death 113 Groups All Years (ICD-10)'
 		DIABETES		= 'Diabetes Flag from Multiple Cause of Death (MCOD)'
		HYPERTEN		= 'Hypertension Flag from Multiple Cause of Death (MCOD)'
		HIPFRACT		= 'Hip Fracture Flag from Multiple Cause of Death (MCOD)'
		WGT_NEW			= "Weight Adjusted for Ineligible Respondents Person-level Sample Weight"
		SA_WGT_NEW		= "Weight Adjusted for Ineligible Respondents Sample Adult Sample Weight"

     ;

   * ASSOCIATE VARIABLES WITH FORMAT VALUES;
   FORMAT    
		ELIGSTAT 		ELIGFMT.          
		MORTSTAT 		MORTFMT.
     	MORTSRCE_NDI 	MRSRCFMT.
     	MORTSRCE_SSA 	MRSRCFMT.
     	MORTSRCE_CMS 	MRSRCFMT.
		CAUSEAVL 		CAUSEFMT.
     	DODQTR   		QRTFMT.           
		DODYEAR  		DODYFMT.
     	DIABETES 		FLAGFMT.          
		HYPERTEN 		FLAGFMT. 
     	HIPFRACT 		FLAGFMT. ;
RUN;