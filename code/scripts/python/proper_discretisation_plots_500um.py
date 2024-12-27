import numpy as np
import  matplotlib.pyplot as plt

dendritic_length = 100

############### mRNA ##################
mRNA_means_100 = np.array([0.0734452,0.0723664,0.0713034,0.070256,0.069224,0.0682072,0.0672053,0.0662181,0.0652454,0.064287,0.0633427,0.0624123,0.0614955,0.0605922,0.0597021,0.0588252,0.0579611,0.0571097,0.0562708,0.0554442,0.0546298,0.0538274,0.0530367,0.0522576,0.05149,0.0507337,0.0499884,0.0492542,0.0485307,0.0478178,0.0471154,0.0464233,0.0457414,0.0450695,0.0444075,0.0437552,0.0431125,0.0424792,0.0418552,0.0412404,0.0406346,0.0400377,0.0394496,0.0388701,0.0382992,0.0377366,0.0371823,0.0366361,0.0360979,0.0355677,0.0350452,0.0345305,0.0340232,0.0335235,0.033031,0.0325459,0.0320678,0.0315967,0.0311326,0.0306753,0.0302247,0.0297807,0.0293433,0.0289123,0.0284876,0.0280691,0.0276568,0.0272506,0.0268503,0.0264559,0.0260673,0.0256844,0.0253071,0.0249353,0.0245691,0.0242082,0.0238526,0.0235022,0.023157,0.0228168,0.0224817,0.0221514,0.021826,0.0215054,0.0211896,0.0208783,0.0205716,0.0202695,0.019972,0.0196794,0.0193928,0.0191165,0.0188652,0.0186892,0.0187618,0.0196773,0.0234751,0.0371539,0.084728,0.2486,0.811521])
n_compartments = len(mRNA_means_100) 
density = mRNA_means_100/dendritic_length*n_compartments

plt.plot(np.linspace(0, dendritic_length, n_compartments), density)


mRNA_means_200 = np.array([0.03706,0.0367853,0.0365126,0.036242,0.0359734,0.0357068,0.0354421,0.0351794,0.0349186,0.0346598,0.0344029,0.0341479,0.0338948,0.0336436,0.0333942,0.0331467,0.032901,0.0326572,0.0324151,0.0321749,0.0319364,0.0316997,0.0314647,0.0312315,0.031,0.0307702,0.0305422,0.0303158,0.0300911,0.029868,0.0296467,0.0294269,0.0292088,0.0289923,0.0287774,0.0285641,0.0283524,0.0281423,0.0279337,0.0277266,0.0275211,0.0273171,0.0271147,0.0269137,0.0267142,0.0265162,0.0263197,0.0261246,0.0259309,0.0257387,0.025548,0.0253586,0.0251706,0.0249841,0.0247989,0.0246151,0.0244326,0.0242515,0.0240718,0.0238934,0.0237163,0.0235405,0.023366,0.0231928,0.0230209,0.0228503,0.0226809,0.0225128,0.0223459,0.0221803,0.0220159,0.0218527,0.0216908,0.02153,0.0213704,0.021212,0.0210548,0.0208987,0.0207438,0.0205901,0.0204375,0.020286,0.0201356,0.0199864,0.0198382,0.0196912,0.0195452,0.0194004,0.0192566,0.0191138,0.0189722,0.0188316,0.018692,0.0185534,0.0184159,0.0182794,0.0181439,0.0180094,0.017876,0.0177435,0.0176119,0.0174814,0.0173518,0.0172232,0.0170956,0.0169689,0.0168431,0.0167182,0.0165943,0.0164713,0.0163492,0.0162281,0.0161078,0.0159884,0.0158699,0.0157523,0.0156355,0.0155196,0.0154046,0.0152904,0.0151771,0.0150646,0.0149529,0.0148421,0.0147321,0.0146229,0.0145145,0.0144069,0.0143001,0.0141941,0.0140889,0.0139845,0.0138809,0.013778,0.0136758,0.0135745,0.0134739,0.013374,0.0132749,0.0131765,0.0130788,0.0129819,0.0128856,0.0127901,0.0126953,0.0126012,0.0125078,0.0124151,0.0123231,0.0122318,0.0121411,0.0120511,0.0119618,0.0118731,0.0117851,0.0116978,0.0116111,0.011525,0.0114396,0.0113548,0.0112706,0.0111871,0.0111042,0.0110219,0.0109402,0.0108591,0.0107786,0.0106987,0.0106194,0.0105407,0.0104626,0.010385,0.0103081,0.0102317,0.0101558,0.0100805,0.0100058,0.00993167,0.00985805,0.00978498,0.00971246,0.00964047,0.00956903,0.00949812,0.00942777,0.00935803,0.00928905,0.00922123,0.00915568,0.00909542,0.00904864,0.00903761,0.00912275,0.00946816,0.0105194,0.013486,0.0216523,0.0439349,0.104543,0.269207,0.716383])
n_compartments = len(mRNA_means_200)
density = mRNA_means_200/dendritic_length*n_compartments

plt.plot(np.linspace(0, dendritic_length, n_compartments), density)


plt.show()

############### Proteins ##################
protein_means_100 = np.array([871.418,871.367,871.268,871.122,870.933,870.702,870.432,870.124,869.781,869.405,868.998,868.562,868.098,867.61,867.099,866.566,866.014,865.444,864.858,864.258,863.646,863.023,862.392,861.753,861.108,860.459,859.808,859.155,858.503,857.853,857.207,856.566,855.931,855.304,854.686,854.079,853.484,852.902,852.334,851.782,851.248,850.732,850.236,849.76,849.307,848.876,848.471,848.091,847.737,847.412,847.115,846.849,846.614,846.411,846.241,846.106,846.006,845.943,845.917,845.93,845.983,846.075,846.21,846.387,846.607,846.871,847.181,847.538,847.941,848.393,848.893,849.444,850.046,850.699,851.405,852.164,852.978,853.848,854.773,855.756,856.797,857.896,859.055,860.275,861.556,862.899,864.305,865.775,867.31,868.911,870.578,872.313,874.115,875.987,877.928,879.938,882.017,884.155,886.324,888.419,890.082])
n_compartments = len(protein_means_100)
density = protein_means_100/dendritic_length*n_compartments

plt.plot(np.linspace(0, dendritic_length, n_compartments), density)


protein_means_200 = np.array([437.847,437.84,437.827,437.808,437.783,437.753,437.716,437.674,437.626,437.573,437.515,437.452,437.384,437.311,437.234,437.152,437.066,436.975,436.881,436.782,436.68,436.574,436.464,436.351,436.235,436.115,435.992,435.867,435.738,435.607,435.473,435.337,435.198,435.057,434.914,434.769,434.622,434.473,434.322,434.17,434.016,433.861,433.705,433.548,433.389,433.23,433.069,432.908,432.746,432.584,432.422,432.258,432.095,431.932,431.768,431.605,431.441,431.278,431.116,430.953,430.791,430.63,430.47,430.31,430.151,429.993,429.836,429.681,429.526,429.373,429.221,429.071,428.922,428.775,428.63,428.486,428.345,428.205,428.068,427.932,427.799,427.668,427.54,427.414,427.29,427.169,427.051,426.936,426.823,426.713,426.607,426.503,426.403,426.305,426.212,426.121,426.034,425.95,425.87,425.794,425.721,425.652,425.587,425.526,425.468,425.415,425.366,425.321,425.281,425.245,425.213,425.185,425.163,425.144,425.131,425.122,425.118,425.118,425.124,425.135,425.15,425.171,425.197,425.228,425.264,425.306,425.353,425.406,425.464,425.528,425.597,425.672,425.753,425.84,425.932,426.031,426.135,426.246,426.362,426.485,426.614,426.75,426.891,427.04,427.194,427.355,427.523,427.697,427.878,428.066,428.261,428.462,428.67,428.885,429.108,429.337,429.574,429.817,430.068,430.326,430.592,430.865,431.145,431.433,431.729,432.032,432.343,432.661,432.988,433.322,433.664,434.014,434.372,434.738,435.112,435.494,435.885,436.284,436.691,437.106,437.53,437.962,438.403,438.853,439.311,439.778,440.253,440.737,441.231,441.733,442.244,442.764,443.293,443.831,444.377,444.932,445.494,446.058,446.612,447.123,447.5])
n_compartments = len(protein_means_200)
density = protein_means_200/dendritic_length*n_compartments

plt.plot(np.linspace(0, dendritic_length, n_compartments), density)

plt.show()

### Correlations ###
protein_PCCs_100 = np.array([1,0.93141,0.931391,0.931352,0.93129,0.931204,0.931092,0.930953,0.930788,0.930595,0.930375,0.930127,0.929853,0.929551,0.929222,0.928867,0.928484,0.928074,0.927637,0.927173,0.926682,0.926164,0.925618,0.925045,0.924444,0.923815,0.923158,0.922472,0.921757,0.921013,0.920238,0.919434,0.918599,0.917733,0.916835,0.915905,0.914942,0.913945,0.912915,0.91185,0.91075,0.909615,0.908443,0.907233,0.905987,0.904701,0.903377,0.902014,0.90061,0.899165,0.897679,0.896151,0.894581,0.892967,0.891309,0.889607,0.88786,0.886068,0.88423,0.882346,0.880414,0.878436,0.87641,0.874337,0.872215,0.870044,0.867825,0.865557,0.863239,0.860872,0.858456,0.85599,0.853474,0.850908,0.848292,0.845627,0.842912,0.840148,0.837335,0.834472,0.831561,0.828601,0.825593,0.822537,0.819434,0.816284,0.813089,0.809848,0.806562,0.803233,0.799861,0.796447,0.792993,0.789499,0.785968,0.782403,0.778808,0.7752,0.771626,0.768247,0.765604])
n_compartments = len(protein_PCCs_100)
plt.plot(np.linspace(0, dendritic_length, n_compartments), protein_PCCs_100, label="100")


protein_PCCs_200 = np.array([1,0.871996,0.871996,0.871995,0.871991,0.871983,0.871972,0.871957,0.871936,0.871912,0.871882,0.871847,0.871807,0.871761,0.87171,0.871653,0.871591,0.871524,0.87145,0.871371,0.871286,0.871196,0.871099,0.870998,0.87089,0.870776,0.870657,0.870532,0.870402,0.870266,0.870124,0.869976,0.869822,0.869663,0.869498,0.869327,0.869151,0.868969,0.868781,0.868587,0.868388,0.868183,0.867971,0.867755,0.867532,0.867303,0.867069,0.866828,0.866582,0.866329,0.866071,0.865807,0.865536,0.86526,0.864977,0.864688,0.864393,0.864092,0.863784,0.86347,0.86315,0.862823,0.862489,0.862149,0.861803,0.86145,0.86109,0.860723,0.860349,0.859969,0.859581,0.859187,0.858785,0.858377,0.857961,0.857537,0.857107,0.856669,0.856223,0.85577,0.855309,0.854841,0.854365,0.853881,0.853389,0.852889,0.852381,0.851864,0.85134,0.850807,0.850266,0.849716,0.849158,0.848591,0.848016,0.847431,0.846838,0.846236,0.845625,0.845005,0.844376,0.843737,0.843089,0.842432,0.841766,0.841089,0.840404,0.839708,0.839003,0.838288,0.837564,0.836829,0.836084,0.835329,0.834564,0.833789,0.833003,0.832207,0.831401,0.830584,0.829757,0.828919,0.828071,0.827212,0.826342,0.825461,0.824569,0.823667,0.822753,0.821829,0.820893,0.819946,0.818988,0.818019,0.817039,0.816047,0.815045,0.81403,0.813005,0.811968,0.810919,0.809859,0.808788,0.807705,0.806611,0.805505,0.804387,0.803258,0.802117,0.800965,0.799801,0.798625,0.797438,0.796239,0.795029,0.793807,0.792573,0.791328,0.790071,0.788802,0.787522,0.786231,0.784927,0.783613,0.782287,0.780949,0.7796,0.77824,0.776868,0.775486,0.774092,0.772686,0.77127,0.769843,0.768405,0.766956,0.765496,0.764025,0.762544,0.761052,0.75955,0.758038,0.756516,0.754983,0.753441,0.751889,0.750327,0.748756,0.747175,0.745586,0.743988,0.742381,0.740766,0.739143,0.737513,0.735878,0.734242,0.732619,0.731041,0.7296,0.72854])
n_compartments = len(protein_PCCs_200)
plt.plot(np.linspace(0, dendritic_length, n_compartments), protein_PCCs_200, label="200")


plt.legend()

plt.show()