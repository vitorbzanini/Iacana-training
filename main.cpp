#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <bits/stdc++.h>
#include "reliability.h"

#define DIM_X 20
#define DIM_Y 20
#define N_STATES 35
#define N_TASKTYPE  3
#define SYSTEM_SIZE (DIM_X*DIM_Y)  
#define THERMAL_NODES (SYSTEM_SIZE*4)+12  // 4 thermal nodes for each PE plus 12 from the environment

struct Tasks{
    int id,type;
    float power;
    int totalTime,taskSlot,fit,temp,currenty_time;
};

Tasks many_core [DIM_Y][DIM_X];

Tasks tasks [16] =  {{0,2,0.350,2000,0},{1,1,0.242,2000,0},{2,0,0.190,2000,0},{3,0,0.172,2000,0},
                    {4,0,0.167,2000,0},{5,0,0.162,2000,0},{6,0,0.154,2000,0},{7,2,0.359,2000,0},
                    {8,2,0.350,2000,0},{9,1,0.242,2000,0},{10,0,0.190,2000,0},{11,0,0.172,2000,0},
                    {12,0,0.167,2000,0},{13,0,0.162,2000,0},{14,0,0.154,2000,0},{15,2,0.359,2000,0}};

struct floorplan_structure floorplan; 
extern struct UnitRel rel_unit[TOTAL_STRUCTURES];

//struct UnitRel rel_unit[TOTAL_STRUCTURES]; /* Struct containing all reliability data*/

double Binv[THERMAL_NODES][SYSTEM_SIZE];
double Cexp[THERMAL_NODES][THERMAL_NODES];

int power[DIM_Y][DIM_X],current_task_allocated=0;
double power_trace[SYSTEM_SIZE];
double t_steady[THERMAL_NODES];

double TempTraceEnd[THERMAL_NODES];
double Tsteady[THERMAL_NODES];
double Tdifference[THERMAL_NODES];
int SystemFIT[DIM_X*DIM_Y];

int getSouth(int x, int y){
    //printf("\nporraadas\n");
    if(y > 0){
        /*if(makeAddress(x,y-1) == GLOBAL_MASTER_ADDR) return 2;
        else*/ return(many_core[y-1][x].type);
    } else {
        return(-1);
    }
}

int getNorth(int x, int y){
    if(y < DIM_Y-1){
        /*if(makeAddress(x,y+1) == GLOBAL_MASTER_ADDR) return 2;
        else*/ return(many_core[y+1][x].type);
    } else {
        return(-1);
    }
}

int getEast(int x, int y){
    if(x < DIM_X-1){
        /*if(makeAddress(x+1,y) == GLOBAL_MASTER_ADDR) return 2;
        else*/ return(many_core[y][x+1].type);
    } else {
        return(-1);
    }
}

int getWest(int x, int y){
    if(x > 0){
        /*if(makeAddress(x-1,y) == GLOBAL_MASTER_ADDR) return 2;
        else*/ return(many_core[y][x-1].type);
    } else {
        return(-1);
    }
}

// A utility function to swap two elements 
void swap(int* a, int* b) { 
    int t = *a; 
    *a = *b; 
    *b = t; 
}

/* This function takes last element as pivot, places 
the pivot element at its correct position in sorted 
array, and places all smaller (smaller than pivot) 
to left of pivot and all greater elements to right 
of pivot */
int partition (int arr[], int arr2[], int low, int high) { 
    int pivot = arr[high]; // pivot 
    int i = (low - 1); // Index of smaller element and indicates the right position of pivot found so far
    int j;
    for (j = low; j <= (high - 1); j++) { 
        // If current element is smaller than the pivot 
        if (arr[j] < pivot) { 
            i++; // increment index of smaller element 
            swap(&arr[i], &arr[j]);
            swap(&arr2[i], &arr2[j]);
        } 
    } 
    swap(&arr[i + 1], &arr[high]); 
    swap(&arr2[i + 1], &arr2[high]); 
    return (i + 1); 
} 

/* The main function that implements QuickSort 
arr[] --> Array to be sorted, 
low --> Starting index, 
high --> Ending index */
void quickSort(int arr[], int arr2[], int low, int high){ 
    if (low < high){ 
        /* pi is partitioning index, arr[p] is now 
        at right place */
        int pi = partition(arr, arr2, low, high); 
  
        // Separately sort elements before 
        // partition and after partition 
        quickSort(arr, arr2, low, (pi - 1)); 
        quickSort(arr, arr2, (pi + 1), high); 
    }
}

//The place of x will be informed by sucessives sums and the y by num(addr%DIM_X or DIM_Y)
unsigned int API_getPEState(unsigned int id, unsigned int excludeAddr){
    unsigned int ex, ey,x = id%DIM_X,y = (int)id/DIM_Y;
    int state_x, state_y, z, state, a;
    //unsigned int xd, yd, zd, state_xd, state_yd, stated; //, b0, c0, ad, bd, cd, i;

    unsigned int immediate[3];
    //unsigned int diagonal[N_TASKTYPE]; //+

    //15 -> [3][3] 15%
    //Decrementing is the easiest way to make that, i think
    //The place of x will be informed by sucessives sums and the y by num(addr%DIM_X or DIM_Y)
    if(excludeAddr != -1){
        //ex = getXpos(excludeAddr);
        //ey = getYpos(excludeAddr);
    } else {
        ex = -1;
        ey = -1;
    }

    for(a = 0; a < 3; a++){
        immediate[a] = 0;
        //diagonal[a] = 0; //+
    }
    // SOUTH
    if(getSouth(x, y) != -1 && !(x == ex && (y-1) == ey ))
        immediate[getSouth(x,y)]++;
    // NORTH
    if(getNorth(x, y) != -1 && !(x == ex && (y+1) == ey ))
        immediate[getNorth(x,y)]++;
    // WEST
    if(getWest(x, y) != -1 && !((x-1) == ex && y == ey ))
        immediate[getWest(x,y)]++;
    // EAST
    if(getEast(x, y) != -1 && !((x+1) == ex && y == ey ))
        immediate[getEast(x,y)]++;
    
    /*
    // SouthWest
    if(getSouthWest(x,y) != -1 && !((x-1) == ex && (y-1) == ey ))
        diagonal[getSouthWest(x,y)]++;
    // SouthEast
    if(getSouthEast(x,y) != -1 && !((x+1) == ex && (y-1) == ey ))
        diagonal[getSouthEast(x,y)]++;
    // NorthWest
    if(getNorthWest(x, y) != -1 && !((x-1) == ex && (y+1) == ey ))
        diagonal[getNorthWest(x,y)]++;
    // NorthEast
    if(getNorthEast(x, y) != -1 && !((x+1) == ex && (y+1) == ey ))
        diagonal[getNorthEast(x,y)]++;*/     
    
    x = immediate[0];
    y = immediate[1];
    z = immediate[2];
    //xd = diagonal[0];
    //yd = diagonal[1];
    //zd = diagonal[2];

    state_x = (int)(x ? ((x*x*x - 18*x*x + 107*x) / 6) : 0);
    state_y = (int)(y ? ((11*y - y*y - 2*x*y) / 2) : 0);
    state = state_x + state_y + z;

    /*state_xd = (int)(xd ? ((xd*xd*xd - 18*xd*xd + 107*xd) / 6) : 0);
    state_yd = (int)(yd ? ((11*yd - yd*yd - 2*xd*yd) / 2) : 0);
    stated = state_xd + state_yd + zd;

    state = stated*35 + state;*/
    
    /*i=0;
    for(a0=0; a0<5; a0++){
        for(b0=0; b0<5; b0++){
            for(c0=0; c0<5; c0++){
                for(ad=0; ad<5; ad++){
                    for(bd=0; bd<5; bd++){
                        for(cd=0; cd<5; cd++){
                            if(a0 == x && b0 == y && c0 == z && ad == xd && bd == yd && cd == zd)
                                return i;
                            else if( (a0 + b0 + c0) <= 4 && (ad + bd + cd) <= 4)
                                i+=1;
                        }
                    }
                }
            }
        }
    }*/

    if(state >= 35) printf("ERRO CALCULANDO ESTADO: ", state);
    return(state);
}

unsigned int API_GetTaskSlotFromTile(unsigned int addr, unsigned int app){
    //printf("\naddr = %d, app = %d, taskSlot = %d\n", addr, app, many_core[(int)addr/DIM_X][addr%DIM_X].taskSlot);
    //getchar();
    if(many_core[(int)addr/DIM_X][addr%DIM_X].taskSlot > 0){
        many_core[(int)addr/DIM_X][addr%DIM_X].taskSlot = many_core[(int)addr/DIM_X][addr%DIM_X].taskSlot - 1;
        //printf("ean");
        if(many_core[(int)addr/DIM_X][addr%DIM_X].type == -1){
            //printf("\nmany_core[%d][%d]\n",(int)addr/DIM_X,addr%DIM_X);
            many_core[(int)addr/DIM_X][addr%DIM_X] = tasks[app];
            current_task_allocated++;
        }
        return 1;
    }else {
        //printf("returning erro5\n");
        return -1;
    }
}

int API_getMaxIdxfromRow(float *policyTable, unsigned int row, unsigned int n_collumns, unsigned int n_rows){
    unsigned int max = 0, i;
    for( i = 0; i < n_collumns; i++){
        if( *(policyTable + row*n_collumns + i) >= max){
            max = i;
        }
    }
    return max;
}

void API_PrintScoreTable(float scoreTable[N_TASKTYPE][N_STATES]){
    int i, j;
    printf("scoreTable = { ");
    for(i = 0; i < N_TASKTYPE; i++){
        printf(" {");
        for(j = 0; j < N_STATES; j++){
            printf("%d",(int)(scoreTable[i][j]*1000));
            if(j != N_STATES-1){
                printf(", ");
            }
        }
        printf(" }");
        if (i != N_TASKTYPE-1){
            printf(",\n");
        }
    }
    printf(" };\n");
}

void load_matrices(){
    FILE *binvpointer;
    binvpointer = fopen("20x20/binv.txt","r");
    FILE *cexppointer;
    cexppointer = fopen("20x20/cexp.txt","r");

    char line[1200000];
    char *number;
    int column, row;

    //printf("atira4");

    for (row = 0; row < THERMAL_NODES; row++){
        fgets(line, sizeof(line), binvpointer);
        number = strtok(line, " ");
        for(column = 0; column < SYSTEM_SIZE; column++){
            Binv[row][column] = strtod(number, NULL);
            //printf("%f ", Binv[row][column]); 
            number = strtok(NULL, " ");      
        }
    }


    for (row = 0; row < THERMAL_NODES; row++){
        fgets(line, sizeof(line), cexppointer);
        number = strtok(line, " ");
        for(column = 0; column < THERMAL_NODES; column++){
            Cexp[row][column] = strtod(number, NULL);
            //printf("%f ", Cexp[row][column]); 
            number = strtok(NULL, " ");      
        }
    }

    fclose(binvpointer);
    fclose(cexppointer);

    for(int i=0;i<THERMAL_NODES;i++){
        TempTraceEnd[i] = 313.15; // Kelvin
    }

    int unitc;
	for (unitc = 0; unitc < TOTAL_STRUCTURES; unitc++){
        sprintf(floorplan.units[unitc].name, "p%d", unitc);
        floorplan.units[unitc].height = 0.000194; // mem 8Kb
        floorplan.units[unitc].width = 0.000194; // mem 8Kb

        init(&floorplan, unitc);  /* Initialize structures*/
        fitinit(unitc);           /* Initialize FITS for each structure*/
    }
}

//many_core[i][j].power = tasks[index%DIM_X][j].power * (1+);
void calcula_temp(){ 
    int index = 0;
    for (int yi = 0; yi < DIM_Y; yi++){
	    for(int xi = 0; xi < DIM_X; xi++){
		    int teste = rand()%26;
            //Ver isso da porcentagem, está errado
            power_trace[index] = (many_core[yi][xi].type != -1) ? many_core[yi][xi].power + ((teste%2==0) ? (many_core[yi][xi].power*((float)teste/100)) : (-1*(many_core[yi][xi].power*((float)teste/100)))) : 0.1; // ATENÇÃO NA ORDEM, OBSERVA OS FORS VARIAM PRIMEIRO O X! E DEPOIS O Y!!!
		    //printf("\npower = %f\n",many_core[yi][xi].power);
            //printf("\nteste1 = %f\n",(many_core[yi][xi].power*((float)teste/100)));
            //printf("\nteste2 = %f\n",(-1*(many_core[yi][xi].power*((float)teste/100))));
            //printf("\npower_trace[%d] = %f\n",index, power_trace[index]);
            //getchar();
            index++;
        }
	}

    int i, j;
    double heatContributionPower;

    for(i = 0; i < THERMAL_NODES; i++){

        heatContributionPower = 0;
        for(j = 0; j < SYSTEM_SIZE; j++){
            heatContributionPower += Binv[i][j]*power_trace[j];
        }
        Tsteady[i] = heatContributionPower + 318.15; // soma com Temperatura Ambiente
    }

    for(int k = 0; k < THERMAL_NODES; k++) Tdifference[k] = TempTraceEnd[k] - Tsteady[k];    

    for(int k = 0; k < THERMAL_NODES; k++){
        double sumExponentials = 0;
        for(j = 0; j < THERMAL_NODES; j++){
            sumExponentials += Cexp[k][j] * Tdifference[j];
        }
        TempTraceEnd[k] = Tsteady[k] + sumExponentials;
    }
}

void calcula_fit(){ 
    for (int structures=0; structures < TOTAL_STRUCTURES; structures++){
        /* Calculate FIT value by feeding in each structures temperature, activity
            * factor, processor supply voltage, and processor frequency. */
        //printf("Info: Temp : %f , Power: %f, strct = %d",TempTraceEnd[structures], power_trace[structures],structures);
        allmodels(TempTraceEnd[structures], power_trace[structures], 1.0, 1.0, structures);
	}

    for(int i = 0; i < DIM_Y*DIM_X; i++) {}//SystemFIT[i] = (int)rel_unit[i].ind_inst*100;

    int m, n, i = 0, avgFit = 0, totalFit = 0;
    for (m = 0; m < DIM_X; m++){
        for (n = 0; n < DIM_Y; n++){
            avgFit = many_core[n][m].fit;
            // printsv("avgFit ", avgFit);
            totalFit = many_core[n][m].fit << 5;
            // printsv("totalFit1 ", totalFit);
            totalFit = totalFit - avgFit;
            // printsv("totalFit2 ", totalFit);
            totalFit = totalFit + SystemFIT[i];
            // printsv("totalFit3 ", totalFit);
            many_core[n][m].fit = totalFit >> 5;
            // printsv("FIT ", many_core[n][m].fit);
            i++;
        }
    }
}

//int num = (rand() % (upper - lower + 1)) + lower;

// {123,3244,11233123,12314} {0,1,2,3} - antes

// {123,3244,12314,11233123} {0,1,3,2}

void manyCorePrint(){ 
    printf("\n{");
    for(int i=0;i<DIM_X;++i){
        for(int j=0;j<DIM_Y;++j){
            printf("%d , ",(many_core[i][j].type != -1) ? 1 : 0);
        }
        printf("\n");
    }
    printf("}");
}

int main(int argc, char *argv[]){

    srand(time(0));

    load_matrices();

    //printf("\natira\n");
    //getchar();

    for(int i=0;i<DIM_X;++i){
        for(int j=0;j<DIM_Y;++j){
            //printf("\natira2\n");
            many_core[j][i].id = 0;
            many_core[j][i].type = -1;
            many_core[j][i].power = 0.1;
            many_core[j][i].totalTime = -1;
            many_core[j][i].taskSlot = 1;
            many_core[j][i].currenty_time = 0;
            many_core[j][i].fit = 0;
            many_core[j][i].temp = 0;
            //printf("\atira3\n");
        }
    }

    //getchar();

    float scoreTable[3][35] ={ {158752.0, 151306.0, 141767.0, 101257.0, 30320.0, 150988.0, 139441.0, 128386.0, 46420.0, 136751.0, 133097.0, 87598.0, 128093.0, 48736.0, 40427.0, 156500.0, 143513.0, 129063.0, 71294.0, 140630.0, 133150.0, 114801.0, 135749.0, 118135.0, 104357.0, 147431.0, 141911.0, 108081.0, 140140.0, 130474.0, 118725.0, 151058.0, 138647.0, 125856.0, 81869.0  },
                               {134928.0, 123406.0, 105842.0, 72358.0, 2556.0, 123309.0, 105027.0, 96300.0, 16705.0, 112223.0, 101299.0, 63861.0, 106332.0, 82031.0, 60822.0, 115340.0, 103802.0, 100406.0, 63490.0, 123720.0, 105601.0, 99110.0, 108154.0, 88274.0, 93688.0, 121426.0, 105191.0, 87878.0, 108250.0, 98963.0, 105731.0, 111918.0, 100972.0, 107116.0, 62499.0  },
                               {155664.0, 145902.0, 124081.0, 126903.0, 8709.0, 140733.0, 125602.0, 112981.0, 33636.0, 111199.0, 105202.0, 41548.0, 104730.0, 23392.0, 0.0, 143522.0, 125388.0, 127966.0, 18118.0, 122550.0, 115151.0, 98306.0, 129990.0, 104634.0, 40216.0, 129079.0, 120941.0, 108770.0, 129943.0, 110096.0, 110038.0, 134345.0, 99074.0, 96193.0, 96179.0  } };

    int sorted_addr[DIM_X*DIM_Y], sorted_score[DIM_X*DIM_Y], addr=0,slot=0,taskType=0, current_fit=0,toprint=0;
    int state_last[DIM_X*DIM_Y],starting_fit[DIM_X*DIM_Y],state_stability[DIM_X*DIM_Y];

    FILE *fl,*fpower,*fp;
    fl = fopen("FITlog.tsv", "w");
    fp = fopen("SystemTemperature.tsv", "w");
    fpower = fopen("SystemPower.tsv", "w");
    fprintf(fp, "time");
    fprintf(fl, "time");
    fprintf(fpower, "time");
    for(int i=0;i<DIM_X*DIM_Y;i++){
        fprintf(fp, "\t%d",i);
        fprintf(fl, "\t%d",i);
        fprintf(fpower, "\t%d",i);
    }
    fprintf(fp, "\n");
    fprintf(fl, "\n");
    fprintf(fpower, "\n");
    fclose(fp);
    fclose(fl);
    fclose(fpower);

    int num_tasks=15,current_xtask=0,current_ytask=0,current_task=0,k=0,cont=0,random_app=0;

    //---------------------------------------
    // --------- Q-learning stuff -----------
    int tp, state, maxid;
    // Hyperparameters
    float epsilon = 0.1;
    float alpha = 0.01;
    float gamma = 0.6;
    float oldvalue, maxval, reward, delta;


    while(1){
        printf("\ncont = %d\n",cont);
        cont++;
        if(cont == 10000) break; 
        
        if(cont>20){
            if(current_task_allocated <= (int)((DIM_X*DIM_Y)*0.7)){ 
                for(int num = (rand() % 6) + 2; num>1;num--){               
                    for(int i = 0; i < DIM_X; i++){
                        for(int j = 0; j < DIM_Y; j++){
                            //Usa a função do makeAdress para divir a matriz em vetor
                            sorted_addr[k] = k;
                            printf("\nsorted_addr[%d] = %d\n",k,sorted_addr[k]);
                            //para o score_table preciso apenas do taskType e State(getState)
                            sorted_score[k] = scoreTable[tasks[num].type][API_getPEState(k, -1)]*1000;
                            printf("\nsorted_score[%d] = %d\n",k,sorted_score[k]);
                            k++;
                        }
                    }

                    k=0;
                    quickSort(sorted_score, sorted_addr, 0, DIM_X*DIM_Y);
                    random_app = (rand() % 16);
                    
                    //printf("randapp = %d", random_app);
                    if((int)(epsilon*100) > random()%100){ // try something new
                        //////////////////////////////////////////////
                        // gets a random tile 
                        do{
                            printf("\nrandom\n");
                            addr = sorted_addr[random()%(DIM_X*DIM_Y)];//pensar em algo com o num_tasks para gerar os endereços 
                            slot = API_GetTaskSlotFromTile(addr, random_app);
                        }while (slot == -1);
                    } else{ // uses the learned information
                        // try to get the best tile slot
                        for(int i = (DIM_X*DIM_Y)-1; i >= 0; i--){
                            addr = sorted_addr[i];
                            printf("\nantes, addr = %d, i = %d\n", addr, i);  
                            //printf("\n explodes2 \n");
                            slot = API_GetTaskSlotFromTile(addr, random_app);
                            //printf("\n explodes3 \n");
                            if (slot != -1) break;
                        }
                    }
                } 
            }


            if(((cont%200)+1) == 200){ 
                
                for(int i = 0; i < DIM_X*DIM_Y; i++){

                    taskType = many_core[(int)i/DIM_X][i%DIM_X].type;
                    if(taskType != -1){ 
                        state = API_getPEState(i, -1);
                        
                        if(state != state_last[i]){ 
                            state_last[i] = state;
                            starting_fit[i] = many_core[(int)i/DIM_X][i%DIM_X].fit;
                            state_stability[i] = 0;
                            current_fit = many_core[(int)i/DIM_X][i%DIM_X].fit;

                        } else if ( state_stability[i] < 50 ) state_stability[i]++;
                        else if ( state == state_last[i] && state_stability[i] >= 50 ) {
                            // calculates the reward
                            delta = (float)((float)(current_fit/100) - (float)(starting_fit[i]/100));
                            if(delta > 0)
                                printf("delta FIT: +", (int)delta);
                            else
                                printf("delta FIT: -", (int)(delta*-1));
                            //reward =  (Q_rsqrt(7000+(delta*delta)) * delta * -500) + 500;
                            reward =  (1/sqrt((10000+(delta*delta)) * delta * -100) + 100);
                            printf("state ", state, "woned reward: ", (int)reward);

                            // gets the old value
                            oldvalue = scoreTable[taskType][state];

                            // gets the max value from the table
                            maxid = API_getMaxIdxfromRow(&(scoreTable[0][0]), taskType, N_STATES, N_TASKTYPE);
                            maxval = scoreTable[taskType][maxid];

                            // updates the score table
                            scoreTable[taskType][state] = (1 - alpha) * oldvalue + alpha * ( reward + gamma * maxval);
                            //scoreTimes[taskType][state] = scoreTimes[taskType][state] + 1;

                            // saves the current FIT for the next update
                            state_stability[i] = 0;
                            starting_fit[i] = many_core[(int)i/DIM_X][i%DIM_X].fit;

                            // print score table
                            toprint++;
                            if(toprint > 100){
                                API_PrintScoreTable(scoreTable);
                                //API_PrintScoreTimes(scoreTimes);
                                toprint=0;
                            }
                        }
                    }
                }
            }


            printf("\nsaiu daqui2\n");

            calcula_temp();
            calcula_fit();        

            for(int i=0;i<DIM_Y*DIM_Y;++i){
                //printf("\nteste\n");
                //printf("\nprintando o FIT: %f\n",rel_unit[i].ind_inst);
                //getchar();
            }

            FILE *powerlog,*fitlog,*fp;
            if(cont%10==0){
                //printf("\nexecutando, cont = %d e %f\n",cont,((float)cont/1000));
                powerlog = fopen("SystemPower.tsv", "a");
                fp = fopen("SystemTemperature.tsv", "a");
                fitlog = fopen("FITlog.tsv", "a");
                fprintf(fitlog, "%.4f", ((float)cont/1000));
                fprintf(powerlog, "%.4f", ((float)cont/1000));
                fprintf(fp, "%.4f", ((float)cont/1000));
            }

            for(int i=0;i<DIM_Y*DIM_Y;++i){
                    if(cont%10==0){
                        fprintf(powerlog,"\t%f",power_trace[i]);
                        fprintf(fp, "\t%.2f", (((float)(TempTraceEnd[i]*100)/100)-273.15));
                        fprintf(fitlog,"\t%f",rel_unit[i].ind_inst);
                        //printf("\nprintando o FIT: %f\n",rel_unit[i].ind_inst);
                    }
                    //printf("many_core[%d][%d] porran\n", i,j);
                    if(many_core[(int)i/DIM_X][i%DIM_X].type != -1) { 
                        many_core[(int)i/DIM_X][i%DIM_X].currenty_time++; 
                    } 
                    if(many_core[(int)i/DIM_X][i%DIM_X].currenty_time >= many_core[(int)i/DIM_X][i%DIM_X].totalTime && many_core[(int)i/DIM_X][i%DIM_X].type != -1) {
                        printf("\n\ncurrent = %d, other = %d\n",current_task_allocated,(int)((DIM_X*DIM_Y)*0.7));
                        manyCorePrint();
                        //getchar();
                        many_core[(int)i/DIM_X][i%DIM_X].id = 0;
                        many_core[(int)i/DIM_X][i%DIM_X].type = -1;
                        many_core[(int)i/DIM_X][i%DIM_X].power = 0.1;
                        many_core[(int)i/DIM_X][i%DIM_X].totalTime = -1;
                        many_core[(int)i/DIM_X][i%DIM_X].taskSlot = 1;
                        many_core[(int)i/DIM_X][i%DIM_X].currenty_time = 0;
                        many_core[(int)i/DIM_X][i%DIM_X].fit = 0;
                        many_core[(int)i/DIM_X][i%DIM_X].temp = 0;
                        many_core[(int)i/DIM_X][i%DIM_X].taskSlot = 1;
                        many_core[(int)i/DIM_X][i%DIM_X].type = -1;
                        many_core[(int)i/DIM_X][i%DIM_X].currenty_time = 0;
                        current_task_allocated--;
                    }
             }
            if(cont%10==0){
                fprintf(powerlog, "\n");
                fclose(powerlog);
                fprintf(fp, "\n");
                fclose(fp);
                fprintf(fitlog, "\n");
                fclose(fitlog);
            }
        }
    }  
    return 0;
}