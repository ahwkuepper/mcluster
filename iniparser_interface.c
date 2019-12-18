
#include <stdio.h>
#include <stdlib.h>
#include "iniparser.h"

void config_getdouble_(double* value, double *defaultValue, const char key[100], const char conf[1000]) {
	dictionary * dict = iniparser_load(conf);

	if (dict == NULL) {
		*value = *defaultValue;
	} else {
		*value = iniparser_getdouble(dict, key, *defaultValue);
		iniparser_freedict(dict);
	}
}

void config_getint_(int* value, int *defaultValue, char key[100], const char conf[1000]) {
	dictionary * dict = iniparser_load(conf);
	if (dict == NULL) {
		*value = *defaultValue;
	} else {
		*value = iniparser_getint(dict, key, *defaultValue);
		iniparser_freedict(dict);
	}
}

void config_getstr_(char value[500], char defaultValue[500], char key[1000], const char conf[1000]) {
	dictionary * dict = iniparser_load(conf);
	if (dict == NULL) {
		strncpy(value, defaultValue, 500);
	} else {
		char * str;
		str = iniparser_getstring(dict, key, "NOT VALID KEY");
   		if (str=="NOT VALID KEY") strncpy(value, defaultValue, 500);
		else strncpy(value, str, 500);
		iniparser_freedict(dict);
	}
}

void char_to_arrayint_(char str[], int *array){
	char *pt;
	int i;
	i=0;
	pt = strtok(str, ", ");

	while (pt != NULL) {
		array[i] = atoi(pt);
		pt = strtok (NULL, ", ");
		i++;
	}
}

void char_to_arraydouble_(char str[], double *array){
	char *pt;
	int i;
	i=0;
	pt = strtok(str,", ");
	while (pt != NULL) {
		array[i] = atof(pt);
		pt = strtok (NULL, ", ");
		i++;
	}
}

void config_validate_(const char conf[10000]) {
	dictionary * dict = iniparser_load(conf);

	int i, k, m, secCount, keysCount, found;
	char* secName;
	char** keys;
	const int expectedCount = 30;
	char* expected[expectedCount];

	if (dict == NULL)
		return;

	expected[0] = "mcluster:n";
	expected[1] = "mcluster:fracb";
	expected[2] = "mcluster:initialmodel";
	expected[3] = "mcluster:w0";
	expected[4] = "mcluster:S";
	expected[5] = "mcluster:fractal";
	expected[6] = "mcluster:qvir";
	expected[7] = "mcluster:mfunc";
	expected[8] = "mcluster:single_mass";
	expected[9] = "mcluster:mlow";
	expected[10] = "mcluster:mup";
	expected[11] = "mcluster:alpha_imf";
	expected[12] = "mcluster:mlim_imf";
	expected[13] = "mcluster:alpha_l3";
	expected[14] = "mcluster:beta_l3";
	expected[15] = "mcluster:mu_l3";
	expected[16] = "mcluster:pairing";
	expected[17] = "mcluster:adis";
	expected[18] = "mcluster:eigen";
	expected[19] = "mcluster:amin";
	expected[20] = "mcluster:amax";
	expected[21] = "mcluster:tf";
	expected[22] = "mcluster:rbar";
	expected[23] = "mcluster:rh_mcl";
	expected[24] = "mcluster:conc_pop";
	expected[25] = "mcluster:potential_energy";
	expected[26] = "mcluster:epoch";
	expected[27] = "mcluster:zini";
	expected[28] = "mcluster:seedmc";
	expected[29] = "mcluster:outputf";


	secCount = iniparser_getnsec(dict);
	for (i = 0; i < secCount; ++i) {
		secName = iniparser_getsecname(dict, i);

		keys = iniparser_getseckeys(dict, secName);
		keysCount = iniparser_getsecnkeys(dict, secName);
		for (k = 0; k < keysCount; ++k) {
//			printf("Checking %s...\n", keys[k]);
			found = 0;
			for (m = 0; m < expectedCount; ++m) {
				if (expected[m] != NULL && strcasecmp(keys[k], expected[m]) == 0) {
					found = 1;
					break;
				}
			}

			if (found == 0) {
				printf("WARNING: key %s from mcluster.ini is not recognized! Check the file mocca.ini for errors.\n", keys[k]);
			}
		}
	}

	iniparser_freedict(dict);
}


