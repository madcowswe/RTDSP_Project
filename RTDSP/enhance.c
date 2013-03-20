/*************************************************************************************
			       DEPARTMENT OF ELECTRICAL AND ELECTRONIC ENGINEERING
					   		     IMPERIAL COLLEGE LONDON 

 				      EE 3.19: Real Time Digital Signal Processing
					       Dr Paul Mitcheson and Daniel Harvey

				        		 PROJECT: Frame Processing

 				            ********* ENHANCE. C **********
							 Shell for speech enhancement 

  		Demonstrates overlap-add frame processing (interrupt driven) on the DSK. 

 *************************************************************************************
 				             By Danny Harvey: 21 July 2006
							 Updated for use on CCS v4 Sept 2010
 ************************************************************************************/
/*
 *	You should modify the code so that a speech enhancement project is built 
 *  on top of this template.
 */
/**************************** Pre-processor statements ******************************/
//  library required when using calloc
#include <stdlib.h>
//  Included so program can make use of DSP/BIOS configuration tool.  
#include "dsp_bios_cfg.h"

/* The file dsk6713.h must be included in every program that uses the BSL.  This 
   example also includes dsk6713_aic23.h because it uses the 
   AIC23 codec module (audio interface). */
#include "dsk6713.h"
#include "dsk6713_aic23.h"

// math library (trig functions)
#include <math.h>

/* Some functions to help with Complex algebra and FFT. */
#include "cmplx.h"      
#include "fft_functions.h"  

// Some functions to help with writing/reading the audio ports when using interrupts.
#include <helper_functions_ISR.h>

#define WINCONST 0.85185			/* 0.46/0.54 for Hamming window */
#define FSAMP 8000.0		/* sample frequency, ensure this matches Config for AIC */
#define FFTLEN 256					/* fft length = frame length 256/8000 = 32 ms*/
#define NFREQ (1+FFTLEN/2)			/* number of frequency bins from a real FFT */
#define OVERSAMP 4					/* oversampling ratio (2 or 4) */  
#define FRAMEINC (FFTLEN/OVERSAMP)	/* Frame increment */
#define CIRCBUF (FFTLEN+FRAMEINC)	/* length of I/O buffers */

#define NumMbuff 4

#define OUTGAIN 16000.0				/* Output gain for DAC */
#define INGAIN  (1.0/16000.0)		/* Input gain for ADC  */
// PI defined here for use in your code 
#define PI 3.141592653589793
#define TFRAME (FRAMEINC/FSAMP)       /* time between calculation of each frame */

#define TMrotate 1
#define MrotateFramecount ((int)(TMrotate/TFRAME))

#define ALPHA 6
#define LAMBDA 0.00
#define TauOP1 0.04

float alpha = ALPHA;
float lambda = LAMBDA;
float kop1 = 0.85; //init in init = exp(-TFRAME/TauOP1);
float nonlinclip = 0.75;
float postgain = 1;


#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define min(a,b) \
({ __typeof__ (a) _a = (a); \
   __typeof__ (b) _b = (b); \
 _a < _b ? _a : _b; })

/******************************* Global declarations ********************************/

/* Audio port configuration settings: these values set registers in the AIC23 audio 
   interface to configure it. See TI doc SLWS106D 3-3 to 3-10 for more info. */
DSK6713_AIC23_Config Config = { \
			 /**********************************************************************/
			 /*   REGISTER	            FUNCTION			      SETTINGS         */ 
			 /**********************************************************************/\
    0x0017,  /* 0 LEFTINVOL  Left line input channel volume  0dB                   */\
    0x0017,  /* 1 RIGHTINVOL Right line input channel volume 0dB                   */\
    0x01f9,  /* 2 LEFTHPVOL  Left channel headphone volume   0dB                   */\
    0x01f9,  /* 3 RIGHTHPVOL Right channel headphone volume  0dB                   */\
    0x0011,  /* 4 ANAPATH    Analog audio path control       DAC on, Mic boost 20dB*/\
    0x0000,  /* 5 DIGPATH    Digital audio path control      All Filters off       */\
    0x0000,  /* 6 DPOWERDOWN Power down control              All Hardware on       */\
    0x0043,  /* 7 DIGIF      Digital audio interface format  16 bit                */\
    0x008d,  /* 8 SAMPLERATE Sample rate control        8 KHZ-ensure matches FSAMP */\
    0x0001   /* 9 DIGACT     Digital interface activation    On                    */\
			 /**********************************************************************/
};

// Codec handle:- a variable used to identify audio interface  
DSK6713_AIC23_CodecHandle H_Codec;

float *inbuffer, *outbuffer;   		/* Input/output circular buffers */
//float *inframe, *outframe;          /* Input and output frames */
float *inwin, *outwin;              /* Input and output windows */
float* powbinstate;
complex* procframe;
complex* procframeprepipe;
float ingain, outgain;				/* ADC and DAC gains */
float cpufrac; 						/* Fraction of CPU time used */
float* powratiobuffs[3];
float* Mbuffs[NumMbuff];
volatile int io_ptr=0;              /* Input/ouput pointer for circular buffers */
volatile int frame_ptr=0;           /* Frame pointer */

 /******************************* Function prototypes *******************************/
void init_hardware(void);    	/* Initialize codec */ 
void init_HWI(void);            /* Initialize hardware interrupts */
void ISR_AIC(void);             /* Interrupt service routine for codec */
void process_frame(void);       /* Frame processing routine */
           
/********************************** Main routine ************************************/
void main()
{      

  	int k, i; // used in various for loops
  
/*  Initialize and zero fill arrays */  

	inbuffer	= (float *) calloc(CIRCBUF, sizeof(float));	/* Input array */
    outbuffer	= (float *) calloc(CIRCBUF, sizeof(float));	/* Output array */
	//inframe		= (float *) calloc(FFTLEN, sizeof(float));	/* Array for processing*/
    //outframe	= (float *) calloc(FFTLEN, sizeof(float));	/* Array for processing*/
    //int(*Mbuff)[FFTLEN] = malloc((sizeof *Mbuff) * NumMbuff);
    procframeprepipe = (complex *) calloc(FFTLEN, sizeof(complex));
    powbinstate	= (float *) calloc(FFTLEN, sizeof(float));
    procframe	= (complex *) calloc(FFTLEN, sizeof(complex));
    inwin		= (float *) calloc(FFTLEN, sizeof(float));	/* Input window */
    outwin		= (float *) calloc(FFTLEN, sizeof(float));	/* Output window */
    
    for (i = 0; i < 3; ++i)
    {
    	powratiobuffs[i] = (float *) calloc(FFTLEN, sizeof(float));
    }

    for (i = 0; i < NumMbuff; ++i)
    {
    	Mbuffs[i] = (float *) calloc(FFTLEN, sizeof(float));
    }
    
    //kop1 = exp(-TFRAME/TauOP1);
	
	/* initialize board and the audio port */
  	init_hardware();
  
  	/* initialize hardware interrupts */
  	init_HWI();    
  
/* initialize algorithm constants */  
                       
  	for (k=0;k<FFTLEN;k++)
	{                           
	inwin[k] = sqrt((1.0-WINCONST*cos(PI*(2*k+1)/FFTLEN))/OVERSAMP);
	outwin[k] = inwin[k]; 
	} 
  	ingain=INGAIN;
  	outgain=OUTGAIN;        

 							
  	/* main loop, wait for interrupt */  
  	while(1) 	process_frame();
}
    
/********************************** init_hardware() *********************************/  
void init_hardware()
{
    // Initialize the board support library, must be called first 
    DSK6713_init();
    
    // Start the AIC23 codec using the settings defined above in config 
    H_Codec = DSK6713_AIC23_openCodec(0, &Config);

	/* Function below sets the number of bits in word used by MSBSP (serial port) for 
	receives from AIC23 (audio port). We are using a 32 bit packet containing two 
	16 bit numbers hence 32BIT is set for  receive */
	MCBSP_FSETS(RCR1, RWDLEN1, 32BIT);	

	/* Configures interrupt to activate on each consecutive available 32 bits 
	from Audio port hence an interrupt is generated for each L & R sample pair */	
	MCBSP_FSETS(SPCR1, RINTM, FRM);

	/* These commands do the same thing as above but applied to data transfers to the 
	audio port */
	MCBSP_FSETS(XCR1, XWDLEN1, 32BIT);	
	MCBSP_FSETS(SPCR1, XINTM, FRM);	
	

}
/********************************** init_HWI() **************************************/ 
void init_HWI(void)
{
	IRQ_globalDisable();			// Globally disables interrupts
	IRQ_nmiEnable();				// Enables the NMI interrupt (used by the debugger)
	IRQ_map(IRQ_EVT_RINT1,4);		// Maps an event to a physical interrupt
	IRQ_enable(IRQ_EVT_RINT1);		// Enables the event
	IRQ_globalEnable();				// Globally enables interrupts

}
        
/******************************** process_frame() ***********************************/  
void process_frame(void)
{
	int k, m, i; 
	int io_ptr0;
	int clearM;

	/* work out fraction of available CPU time used by algorithm */    
	cpufrac = ((float) (io_ptr & (FRAMEINC - 1)))/FRAMEINC;  
		
	/* wait until io_ptr is at the start of the current frame */ 	
	while((io_ptr/FRAMEINC) != frame_ptr); 
	
	/* then increment the framecount (wrapping if required) */ 
	if (++frame_ptr >= (CIRCBUF/FRAMEINC)) frame_ptr=0;
 	
 	/* save a pointer to the position in the I/O buffers (inbuffer/outbuffer) where the 
 	data should be read (inbuffer) and saved (outbuffer) for the purpose of processing */
 	io_ptr0=frame_ptr * FRAMEINC;
	
	/* copy input data from inbuffer into inframe (starting from the pointer position) */ 
	 
	m=io_ptr0;
    for (k=0;k<FFTLEN;k++)
	{                           
		procframe[k] = cmplx(inbuffer[m] * inwin[k], 0);
		if (++m >= CIRCBUF) m=0; /* wrap if required */
	} 
	
	/************************* DO PROCESSING OF FRAME  HERE **************************/
	
	
	/* please add your code, at the moment the code simply copies the input to the 
	ouptut with no processing */

	
	fft(FFTLEN, procframe);
				
    static int MRotate_ctr = 0;

	if (MRotate_ctr == MrotateFramecount){

		MRotate_ctr = 0;

		//rotate Mbuffs
		float* rottemp = Mbuffs[NumMbuff-1];
		for (i = NumMbuff-1; i != 0; --i)
			Mbuffs[i] = Mbuffs[i-1];
		Mbuffs[0] = rottemp;

		//Init new mbuff
		clearM = 1;
	} else {
		MRotate_ctr++;
		clearM = 0;
	}

	powratiobuffs[2] = powratiobuffs[1];
	powratiobuffs[1] = powratiobuffs[0];

	static float lastallbinpower = 0;
	static float lastallbinsigpower = 0;
	float allbinpower = lastallbinpower;
	float allbinsigpower = lastallbinsigpower;
	lastallbinpower = 0;
	lastallbinsigpower = 0;

	for (k = 0; k < FFTLEN; ++k)
	{
		float curramp = cabs(procframe[k]);
		float currpow = curramp * curramp;

		//LPF currpow
		currpow = (1-kop1)*currpow + kop1*powbinstate[k];
		powbinstate[k] = currpow;

		float currnoisebin = clearM ? currpow : min(Mbuffs[0][k], currpow);
		Mbuffs[0][k] = currnoisebin;

		for (i = 1; i < NumMbuff; ++i)
			currnoisebin = min(currnoisebin, Mbuffs[i][k]);

		currnoisebin *= alpha;
		
		lastallbinpower += currnoisebin;
		lastallbinsigpower += currpow;
		
		powratiobuffs[0][k] = currnoisebin/(curramp*curramp);
		//float g = (powratio > (nonlinclip * (allbinpower/**allbinpower*/))) ? 0 : max(lambda, 1-sqrt(powratio));
		//procframe[k] = rmul(g, procframe[k]);
	}

	//swap new and old
	complex* tempframe = procframeprepipe;
	procframeprepipe = procframe;
 	procframe = tempframe;

 	for (k = 0; k < FFTLEN; ++k)
 	{
 		float powratiomin = min(powratiobuffs[0][k], powratiobuffs[2][k]);//, powratiobuffs[1][k]); //NSR!!!

 		float g = (powratiomin > (nonlinclip * (allbinpower/allbinsigpower))) ? 0 : max(lambda, 1-sqrt(powratiobuffs[1][k]));
 		procframe[k] = rmul(g*postgain, procframe[k]);
 	}
 	
	ifft(FFTLEN, procframe);
	
	/********************************************************************************/
	
    /* multiply outframe by output window and overlap-add into output buffer */  
                           
	m=io_ptr0;
    
    for (k=0;k<(FFTLEN-FRAMEINC);k++) 
	{    										/* this loop adds into outbuffer */                       
	  	outbuffer[m] = outbuffer[m] + procframe[k].r *outwin[k];   
		if (++m >= CIRCBUF) m=0; /* wrap if required */
	}         
    for (;k<FFTLEN;k++) 
	{                           
		outbuffer[m] = procframe[k].r * outwin[k];   /* this loop over-writes outbuffer */        
	    m++;
	}	                                   
}        
/*************************** INTERRUPT SERVICE ROUTINE  *****************************/

// Map this to the appropriate interrupt in the CDB file
   
void ISR_AIC(void)
{       
	short sample;
	/* Read and write the ADC and DAC using inbuffer and outbuffer */
	
	sample = mono_read_16Bit();
	inbuffer[io_ptr] = ((float)sample)*ingain;
		/* write new output data */
	mono_write_16Bit((int)(outbuffer[io_ptr]*outgain)); 
	
	/* update io_ptr and check for buffer wraparound */    
	
	if (++io_ptr >= CIRCBUF) io_ptr=0;
}

/************************************************************************************/
