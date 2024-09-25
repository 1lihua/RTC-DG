
#include"read_input.h"

bool Read_inputfile(char* inputfile){
	FILE* fp = fopen(inputfile, "r");
	if (fp == NULL)
	{
		//LOG("Failed to open file: %s", inputfile);
		return false;
	}
	else
	{
		//LOG(0, "Readingparameters from input file \"%s\"\n\n", inputfile);
	}

	// header line
	//	fgets(buf, BUFSIZ, fp);
	//	if (procid == 0)
	//		fputs(buf,stdout);
	char buf[512];
	while (!feof(fp))
	{

		fgets(buf, BUFSIZ, fp);

		// if find two successive # in this line, move to the next line
		if (strstr(buf, "##") != NULL)
		{
			//			LOG(0, "COMMENT: %.30s\n", buf);
		}
		else
		{
			if (strstr(buf, "SI_meter") != NULL) {
				fgets(buf, BUFSIZ, fp);
				sscanf(buf, "%lf", &SI_meter);
			}

			if (strstr(buf, "Frequence range") != NULL) // Initialization instruction line
			{
				fgets(buf, BUFSIZ, fp);
				sscanf(buf, "%lf %lf %lf", &fre_start, &fre_step, &fre_end);


				//sscanf(buf, "%d", &initEM);
			}


			if (strstr(buf, "EM Initialization") != NULL) // Initialization instruction line
			{
				fgets(buf, BUFSIZ, fp);
				//sscanf(buf, "%d", &initEM);
			}
			if (strstr(buf, "LUMP") != NULL) {
				int num_lump=0;
				fgets(buf, BUFSIZ, fp);
				sscanf(buf, "%d", &num_lump);
				H_lump = (double*)calloc(num_lump, sizeof(double));
				W_lump = (double*)calloc(num_lump, sizeof(double));
				for (int n = 0; n < num_lump; n++)
				{
					fgets(buf, BUFSIZ, fp);
					//sscanf(buf, "%lf %lf %lf %lf %lf %lf %lf %lf",  &DIEeps[n], &DIEmu[n], &DIEsigma[n]);
					sscanf(buf, "%lf %lf ", &H_lump[n], &W_lump[n]);
					//LOG(0, "Die Material[%d]: eps %f, mu %f, sigma %f  \n", n, DIEeps[n], DIEmu[n], DIEsigma[n]);
				}
			}
			if (strstr(buf, "WAVE") != NULL) {
				//num_wave;

				fgets(buf, BUFSIZ, fp);
				sscanf(buf, "%lf", &wave_zzz);

				fgets(buf, BUFSIZ, fp);
				sscanf(buf, "%d", &num_wave);

				fgets(buf, BUFSIZ, fp);
				sscanf(buf, "%lf %lf ", &rb_wave, &ra_wave);

				x_wave = (double*)calloc(num_wave, sizeof(double));
				y_wave = (double*)calloc(num_wave, sizeof(double));
				z_wave = (double*)calloc(num_wave, sizeof(double));
				phase_wave = (double*)calloc(num_wave, sizeof(double));
				Amplitude_wave = (double*)calloc(num_wave, sizeof(double));
				for (int n = 0; n < num_wave; n++)
				{
					fgets(buf, BUFSIZ, fp);
					//sscanf(buf, "%lf %lf %lf %lf %lf %lf %lf %lf",  &DIEeps[n], &DIEmu[n], &DIEsigma[n]);
					sscanf(buf, "%lf %lf %lf %lf %lf ", &x_wave[n], &y_wave[n], & z_wave[n], &phase_wave[n], &Amplitude_wave[n]);
					//LOG(0, "Die Material[%d]: eps %f, mu %f, sigma %f  \n", n, DIEeps[n], DIEmu[n], DIEsigma[n]);
				}
			}

			if (strstr(buf, "Materials") != NULL) // Material parameters instruction line
			{
				int MatTypeNum;
				fgets(buf, BUFSIZ, fp);
				sscanf(buf, "%d", &MatTypeNum);
				//LOG(0, "Number of material types: %d \n", MatTypeNum);

				int tread = 0;
				for (int t = 0; t < MatTypeNum; t++)
				{
					fgets(buf, BUFSIZ, fp);
					if (strstr(buf, "Die") != NULL) // Dielectric material, we know there will be a space here
					{
						int DIEnum;
						fgets(buf, BUFSIZ, fp);
						sscanf(buf, "%d", &DIEnum);
						if (DIEnum <= 0 || DIEnum > 100)
						{
							//LOG("!!!!!!!!!! number of materials in each type should be in [1,100] !!!!!!!!!! \n");
							return false;
						}

						// let's predefine 4 values for each type, in InitMat() we may access them
						DIEeps = (double*)calloc(DIEnum, sizeof(double));
						DIEmu = (double*)calloc(DIEnum, sizeof(double));
						DIEsigma = (double*)calloc(DIEnum, sizeof(double));
						DIEk1 = (double*)calloc(DIEnum, sizeof(double));
						DIEk2 = (double*)calloc(DIEnum, sizeof(double));
						DIEk3 = (double*)calloc(DIEnum, sizeof(double));
						DIERho = (double*)calloc(DIEnum, sizeof(double));
						DIESHC = (double*)calloc(DIEnum, sizeof(double));
						DIEQ = (double*)calloc(DIEnum, sizeof(double));
						DIEE = (double*)calloc(DIEnum, sizeof(double));
						DIEnu = (double*)calloc(DIEnum, sizeof(double));
						DIEalpha = (double*)calloc(DIEnum, sizeof(double));

						for (int n = 0; n < DIEnum; n++)
						{
							fgets(buf, BUFSIZ, fp);
							//sscanf(buf, "%lf %lf %lf %lf %lf %lf %lf %lf",  &DIEeps[n], &DIEmu[n], &DIEsigma[n]);
							sscanf(buf, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &DIEeps[n], &DIEmu[n], &DIEsigma[n], &DIEk1[n], &DIEk2[n], &DIEk3[n], &DIERho[n], &DIESHC[n], &DIEQ[n], &DIEE[n], &DIEnu[n], &DIEalpha[n]);
							//LOG(0, "Die Material[%d]: eps %f, mu %f, sigma %f  \n", n, DIEeps[n], DIEmu[n], DIEsigma[n]);
						}
						tread++;
					}
				}
			}

			if (strstr(buf, "PATTERNF") != NULL) // Initialization instruction line
			{
				fgets(buf, BUFSIZ, fp);
				sscanf(buf, "%lf", &fre_pattern);

			}
			if (strstr(buf, "ET") != NULL) // Initialization instruction line
			{
				E_T = true;

				fgets(buf, BUFSIZ, fp);
				sscanf(buf, "%lf", &alepa_T);
			}

			if (strstr(buf, "NEUMANN1") != NULL) // Initialization instruction line
			{
				fgets(buf, BUFSIZ, fp);
				sscanf(buf, "%lf %lf", &Hcvt, &Td);
			}
			if (strstr(buf, "EM_solve") != NULL) // Initialization instruction line
			{
				fgets(buf, BUFSIZ, fp);
				sscanf(buf, "%d ", &EM_solve);
				fgets(buf, BUFSIZ, fp);
				sscanf(buf, "%lf", &fre_E_T_F);
			}
			if (strstr(buf, "TF_only") != NULL) // Initialization instruction line
			{
				TF_only = true;



			}
			if (strstr(buf, "POWER") != NULL) // Initialization instruction line
			{
				fgets(buf, BUFSIZ, fp);
				sscanf(buf, "%lf ", &powerful);
			}



		} // if ##

	} // feof
	fclose(fp);

}