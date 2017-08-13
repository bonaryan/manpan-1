rem model runs
rem  n   diam    d/t   lambda   s/d
rem [1] [2 3 4] [3 5] [1 2 3] [3 4 5]

(
abaqus cae noGUI=Polygoner_ABQS.py -- 1 2 3 1 3
abaqus cae noGUI=Polygoner_ABQS.py -- 1 2 3 1 4
abaqus cae noGUI=Polygoner_ABQS.py -- 1 2 3 1 5
abaqus cae noGUI=Polygoner_ABQS.py -- 1 2 3 2 3
abaqus cae noGUI=Polygoner_ABQS.py -- 1 2 3 2 4
abaqus cae noGUI=Polygoner_ABQS.py -- 1 2 3 2 5
abaqus cae noGUI=Polygoner_ABQS.py -- 1 2 3 3 3
abaqus cae noGUI=Polygoner_ABQS.py -- 1 2 3 3 4
abaqus cae noGUI=Polygoner_ABQS.py -- 1 2 3 3 5
abaqus cae noGUI=Polygoner_ABQS.py -- 1 2 5 1 3
abaqus cae noGUI=Polygoner_ABQS.py -- 1 2 5 1 4
abaqus cae noGUI=Polygoner_ABQS.py -- 1 2 5 1 5
abaqus cae noGUI=Polygoner_ABQS.py -- 1 2 5 2 3
abaqus cae noGUI=Polygoner_ABQS.py -- 1 2 5 2 4
abaqus cae noGUI=Polygoner_ABQS.py -- 1 2 5 2 5
abaqus cae noGUI=Polygoner_ABQS.py -- 1 2 5 3 3
abaqus cae noGUI=Polygoner_ABQS.py -- 1 2 5 3 4
abaqus cae noGUI=Polygoner_ABQS.py -- 1 2 5 3 5
abaqus cae noGUI=Polygoner_ABQS.py -- 1 3 3 1 3
abaqus cae noGUI=Polygoner_ABQS.py -- 1 3 3 1 4
abaqus cae noGUI=Polygoner_ABQS.py -- 1 3 3 1 5
abaqus cae noGUI=Polygoner_ABQS.py -- 1 3 3 2 3
abaqus cae noGUI=Polygoner_ABQS.py -- 1 3 3 2 4
abaqus cae noGUI=Polygoner_ABQS.py -- 1 3 3 2 5
abaqus cae noGUI=Polygoner_ABQS.py -- 1 3 3 3 3
abaqus cae noGUI=Polygoner_ABQS.py -- 1 3 3 3 4
abaqus cae noGUI=Polygoner_ABQS.py -- 1 3 3 3 5
abaqus cae noGUI=Polygoner_ABQS.py -- 1 3 5 1 3
abaqus cae noGUI=Polygoner_ABQS.py -- 1 3 5 1 4
abaqus cae noGUI=Polygoner_ABQS.py -- 1 3 5 1 5
abaqus cae noGUI=Polygoner_ABQS.py -- 1 3 5 2 3
abaqus cae noGUI=Polygoner_ABQS.py -- 1 3 5 2 4
abaqus cae noGUI=Polygoner_ABQS.py -- 1 3 5 2 5
abaqus cae noGUI=Polygoner_ABQS.py -- 1 3 5 3 3
abaqus cae noGUI=Polygoner_ABQS.py -- 1 3 5 3 4
abaqus cae noGUI=Polygoner_ABQS.py -- 1 3 5 3 5
abaqus cae noGUI=Polygoner_ABQS.py -- 1 4 3 1 3
abaqus cae noGUI=Polygoner_ABQS.py -- 1 4 3 1 4
abaqus cae noGUI=Polygoner_ABQS.py -- 1 4 3 1 5
abaqus cae noGUI=Polygoner_ABQS.py -- 1 4 3 2 3
abaqus cae noGUI=Polygoner_ABQS.py -- 1 4 3 2 4
abaqus cae noGUI=Polygoner_ABQS.py -- 1 4 3 2 5
abaqus cae noGUI=Polygoner_ABQS.py -- 1 4 3 3 3
abaqus cae noGUI=Polygoner_ABQS.py -- 1 4 3 3 4
abaqus cae noGUI=Polygoner_ABQS.py -- 1 4 3 3 5
abaqus cae noGUI=Polygoner_ABQS.py -- 1 4 5 1 3
abaqus cae noGUI=Polygoner_ABQS.py -- 1 4 5 1 4
abaqus cae noGUI=Polygoner_ABQS.py -- 1 4 5 1 5
abaqus cae noGUI=Polygoner_ABQS.py -- 1 4 5 2 3
abaqus cae noGUI=Polygoner_ABQS.py -- 1 4 5 2 4
abaqus cae noGUI=Polygoner_ABQS.py -- 1 4 5 2 5
abaqus cae noGUI=Polygoner_ABQS.py -- 1 4 5 3 3
abaqus cae noGUI=Polygoner_ABQS.py -- 1 4 5 3 4
abaqus cae noGUI=Polygoner_ABQS.py -- 1 4 5 3 5
)

rem d	ed/t	lambda	s/d	maxL
rem 500	087.8	0.65	3
rem 700	087.8	0.65	3
rem 900	087.8	0.65	3
rem 500	105.6	0.65	3
rem 700	105.6	0.65	3
rem 900	105.6	0.65	3
rem 500	087.8	1.00	3
rem 700	087.8	1.00	3
rem 900	087.8	1.00	3
rem 500	105.6	1.00	3
rem 700	105.6	1.00	3
rem 900	105.6	1.00	3
rem 500	087.8	1.25	3
rem 700	087.8	1.25	3
rem 900	087.8	1.25	3
rem 500	105.6	1.25	3
rem 700	105.6	1.25	3
rem 900	105.6	1.25	3
rem 500	087.8	0.65	4
rem 700	087.8	0.65	4
rem 900	087.8	0.65	4
rem 500	105.6	0.65	4
rem 700	105.6	0.65	4
rem 900	105.6	0.65	4
rem 500	087.8	1.00	4
rem 700	087.8	1.00	4
rem 900	087.8	1.00	4
rem 500	105.6	1.00	4
rem 700	105.6	1.00	4
rem 900	105.6	1.00	4
rem 500	087.8	1.25	4
rem 700	087.8	1.25	4
rem 900	087.8	1.25	4
rem 500	105.6	1.25	4
rem 700	105.6	1.25	4
rem 900	105.6	1.25	4
rem 500	087.8	0.65	5
rem 700	087.8	0.65	5
rem 900	087.8	0.65	5
rem 500	105.6	0.65	5
rem 700	105.6	0.65	5
rem 900	105.6	0.65	5
rem 500	087.8	1.00	5
rem 700	087.8	1.00	5
rem 900	087.8	1.00	5
rem 500	105.6	1.00	5
rem 700	105.6	1.00	5
rem 900	105.6	1.00	5
rem 500	087.8	1.25	5
rem 700	087.8	1.25	5
rem 900	087.8	1.25	5
rem 500	105.6	1.25	5
rem 700	105.6	1.25	5
rem 900	105.6	1.25	5