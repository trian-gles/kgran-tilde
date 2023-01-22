/**
	@file
	kgran~ - a granulation algorithm designed by Mara Helmuth, rewritten and ported to Max by Kieran McAuliffe

*/

/* 
	should have an inlet optionally used for the phase of the buffer.  If it is not connected, this will be 0
*/

#include "ext.h"
#include "z_dsp.h"
#include "math.h"
#include "ext_buffer.h"
#include "ext_atomic.h"
#include "ext_obex.h"
#include <stdlib.h>
#include <stdbool.h>

#define MAXGRAINS 1000
#define MIDC_OFFSET (261.62556530059868 / 256.0)
#define M_LN2	0.69314718055994529



typedef struct Grain {
	float waveSampInc; 
	float ampSampInc; 
	float ampPhase; 
	float endTime; 
	float panR; 
	float panL; 
	float currTime; 
	bool isplaying;
	} Grain;

typedef struct _kgran {
	t_pxobject w_obj;
	t_buffer_ref *w_buf;
	t_buffer_ref *w_env;
	t_symbol *w_name;
	t_symbol *w_envname;
	
	t_bool running;
	t_bool w_buffer_modified;
	Grain grains[MAXGRAINS];
	
	long w_len;
	long w_envlen;
	
	double transLow;
	double transMid;
	double transHigh;
	double transTight;

	double grainDurLow;
	double grainDurMid;
	double grainDurHigh;
	double grainDurTight;
	
	double grainHeadLow;
	double grainHeadMid;
	double grainHeadHigh;
	double grainHeadTight;

	double panLow;
	double panMid;
	double panHigh;
	double panTight;
	
	double grainRateVarLow;
	double grainRateVarMid;
	double grainRateVarHigh;
	double grainRateVarTight;
	
	int newGrainCounter;
	float grainRate;
	
	double oneover_cpsoct10;
	
	short w_connected;
} t_kgran;




void *kgran_new(t_symbol *s,  long argc, t_atom *argv);
void kgran_free(t_kgran *x);
t_max_err kgran_notify(t_kgran *x, t_symbol *s, t_symbol *msg, void *sender, void *data);
void kgran_assist(t_kgran *x, void *b, long m, long a, char *s);
void kgran_start(t_kgran *x);
void kgran_stop(t_kgran *x);
void kgran_perform64(t_kgran *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam);
void kgran_dsp64(t_kgran *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags);
void kgran_setbuffers(t_kgran* x, t_symbol* s, long ac, t_atom* av);
void kgran_set(t_kgran* x, t_symbol* s, long argc, t_atom* argv);
void kgran_grainrate(t_kgran *x, double rl, double rm, double rh, double rt);
void kgran_graindur(t_kgran *x, double dl, double dm, double dh, double dt);
void kgran_grainhead(t_kgran *x, double hl, double hm, double hh, double ht);
void kgran_trans(t_kgran *x, double fl, double fm, double fh, double ft);
void kgran_pan(t_kgran *x, double pl, double pm, double ph, double pt); 

double rrand() 
{
	double min = -1;
	double max = 1;
    double range = (max - min); 
    double div = RAND_MAX / range;
    return min + (rand() / div);
}

// Taken from RTCMIX code 
// https://github.com/RTcmix/RTcmix/blob/1b04fd3f121a1c65743fde8ea37eb5d65f2cf35c/genlib/pitchconv.c
double octcps(double cps)
{
	return log(cps / MIDC_OFFSET) / M_LN2;
}

double cpsoct(double oct)
{
	return pow(2.0, oct) * MIDC_OFFSET;
}

float oscili(float amp, float si, float *farray, int len, float *phs)
{
	register int i =  *phs;        
	register int k =  (i + 1) % len;  
	float frac = *phs  - i;      
	*phs += si;                 
	while(*phs >= len)
		*phs -= len;       
	return((*(farray+i) + (*(farray+k) - *(farray+i)) *
					   frac) * amp);
}



// From Mara Helmuth
double prob(double low,double mid,double high,double tight)
        // Returns a value within a range close to a preferred value
                    // tightness: 0 max away from mid
                     //               1 even distribution
                      //              2+amount closeness to mid
                      //              no negative allowed
{
	double range, num, sign;

	range = (high-mid) > (mid-low) ? high-mid : mid-low;
	do {
	  	if (rrand() > 0.)
			sign = 1.;
		else  sign = -1.;
	  	num = mid + sign*(pow((rrand()+1.)*.5,tight)*range);
	} while(num < low || num > high);
	return(num);
}

static t_symbol *ps_buffer_modified;
static t_class *s_kgran_class;


void ext_main(void *r)
{
	t_class *c = class_new("kgran~", (method)kgran_new, (method)kgran_free, sizeof(t_kgran), NULL, A_GIMME, 0);

	class_addmethod(c, (method)kgran_dsp64,		"dsp64",	A_CANT, 0);
	class_addmethod(c, (method)kgran_start,		"start", 0);
	class_addmethod(c, (method)kgran_stop,		"stop", 0);
	
	class_addmethod(c, (method)kgran_notify,		"notify",	A_CANT, 0);
	class_addmethod(c, (method)kgran_set, "set", A_GIMME, 0);
	class_addmethod(c, (method)kgran_assist,		"assist",	A_CANT, 0);
	
	// these float methods should be replaced with A_GIMME, see docs
	class_addmethod(c, (method)kgran_grainrate, "grainrate", 
	A_FLOAT, A_FLOAT, A_FLOAT, A_FLOAT, 0);
	
	class_addmethod(c, (method)kgran_graindur, "graindur", 
	A_FLOAT, A_FLOAT, A_FLOAT, A_FLOAT, 0);
	
	class_addmethod(c, (method)kgran_grainhead, "grainhead", 
	A_FLOAT, A_FLOAT, A_FLOAT, A_FLOAT, 0);
	
	class_addmethod(c, (method)kgran_trans, "trans", 
	A_FLOAT, A_FLOAT, A_FLOAT, A_FLOAT, 0);
	
	class_addmethod(c, (method)kgran_pan, "pan", 
	A_FLOAT, A_FLOAT, A_FLOAT, A_FLOAT, 0);

	class_dspinit(c);
	class_register(CLASS_BOX, c);
	s_kgran_class = c;

	ps_buffer_modified = gensym("buffer_modified");
}
/*
	Inlets:
	0 : On/off, grainrate, graindur, trans, pan
	
*/

/* Args:
		p0: wavetable
		p1: grainEnv
	*/
	
// will eventually need to handle for buffers with more than one channel
void *kgran_new(t_symbol *s,  long argc, t_atom *argv)
{
	t_kgran *x = (t_kgran *)object_alloc(s_kgran_class);
	t_symbol *buf=0;
	t_symbol *env=0;

	dsp_setup((t_pxobject *)x,1); // inlet for the buffer head
	buf = atom_getsymarg(0,argc,argv);
	env = atom_getsymarg(1,argc,argv);

	x->w_name = buf;
	x->w_envname = env;
	x->w_buffer_modified = false;
	
	
	//outlets
	outlet_new((t_object *)x, "signal");		// audio outlet l
	outlet_new((t_object *)x, "signal");		// audio outlet r
	
	
	// Setup Grains
	for (size_t i = 0; i < MAXGRAINS; i++){
        x->grains[i] = (Grain){.waveSampInc=0, 
        	.ampSampInc=0,
        	.endTime=0, 
        	.ampPhase=0,
        	.panR=0, 
        	.panL=0, 
        	.currTime=0, 
        	.isplaying=false };
    }

	// create a new buffer reference, initially referencing a buffer with the provided name
	x->w_buf = buffer_ref_new((t_object *)x, x->w_name);
	x->w_env = buffer_ref_new((t_object *)x, x->w_envname);

	t_buffer_obj* b = buffer_ref_getobject(x->w_buf);
	t_buffer_obj* e = buffer_ref_getobject(x->w_env);
	
	x->w_len = buffer_getframecount(b);
	x->w_envlen = buffer_getframecount(e);
	
	x->oneover_cpsoct10 = 1.0 / cpsoct(10.0);

	return (x);
}




void kgran_free(t_kgran *x)
{
	dsp_free((t_pxobject *)x);

	// must free our buffer reference when we will no longer use it
	object_free(x->w_buf);
	object_free(x->w_env);
}


// A notify method is required for our buffer reference
// This handles notifications when the buffer appears, disappears, or is modified.
t_max_err kgran_notify(t_kgran *x, t_symbol *s, t_symbol *msg, void *sender, void *data)
{
	if (msg == ps_buffer_modified) {
		// post("BUFFER CHANGED");
		x->w_buffer_modified = true;
		if (!x->running)
			defer((t_object*)x, (method)kgran_setbuffers, NULL, 0, NULL);
		if (s == x->w_name) {
			return buffer_ref_notify(x->w_buf, s, msg, sender, data);
		}
		else if (s == x->w_envname)
		{
			return buffer_ref_notify(x->w_env, s, msg, sender, data);
		}


	}
		
	
}


////
// SET BUFFER
////

void kgran_setbuffers(t_kgran* x, t_symbol* s, long ac, t_atom* av) {
	object_free(x->w_buf);
	object_free(x->w_env);

	x->w_buf = buffer_ref_new((t_object*)x, x->w_name);
	x->w_env = buffer_ref_new((t_object*)x, x->w_envname);

	t_buffer_obj* b = buffer_ref_getobject(x->w_buf);
	t_buffer_obj* e = buffer_ref_getobject(x->w_env);

	x->w_len = buffer_getframecount(b);
	x->w_envlen = buffer_getframecount(e);
}

void kgran_set(t_kgran* x, t_symbol* s, long argc, t_atom* argv) {
	x->w_name = atom_getsymarg(0, argc, argv);
	x->w_envname = atom_getsymarg(1, argc, argv);
	defer((t_object*)x, (method)kgran_setbuffers, NULL, 0, NULL);
}

void kgran_assist(t_kgran *x, void *b, long m, long a, char *s)
{
	if (m == ASSIST_INLET) {	// inlets
		switch (a) {
		case 0:	snprintf_zero(s, 256, "Various messages");	break;
		}
	}
	else if (m == ASSIST_OUTLET){	// outlet
		switch (a) {
		case 0:	snprintf_zero(s, 256, "(signal) right output");	break;
		case 1:	snprintf_zero(s, 256, "(signal) left output");	break;
		}
	}
}


////
// START AND STOP MSGS
////
void kgran_start(t_kgran *x){
	if (!buffer_ref_exists(x->w_buf) || !buffer_ref_exists(x->w_env))
	{
		error("Make sure you've configured a wavetable buffer and envelope buffer!");
		defer((t_object*)x, (method)kgran_setbuffers, NULL, 0, NULL);
	}

	else
		x->running = true;
}

void kgran_stop(t_kgran *x){
	x->running = false;
}


////
// PARAMETER MESSAGES
////

void kgran_grainrate(t_kgran* x, double rl, double rm, double rh, double rt) {
	x->grainRateVarLow = rl;
	x->grainRateVarMid = fmax(rm, rl);
	x->grainRateVarHigh = fmax(rh, rm);
	x->grainRateVarTight = rt;
}

void kgran_graindur(t_kgran* x, double dl, double dm, double dh, double dt) {
	x->grainDurLow = dl;
	x->grainDurMid = fmax(dm, dl);
	x->grainDurHigh = fmax(dh, dm);
	x->grainDurTight = dt;
}

void kgran_grainhead(t_kgran *x, double hl, double hm, double hh, double ht){
	x->grainHeadLow = fmax(0, hl);
	x->grainHeadMid = fmax(hm, hl);
	x->grainHeadHigh = fmin(fmax(hh, hm), 1);
	x->grainHeadTight = ht;
}

void kgran_trans(t_kgran *x, double fl, double fm, double fh, double ft){
	x->transLow = fl;
	x->transMid = fmax(fm, fl);
	x->transHigh = fmax(fh, fm);
	x->transTight = ft;
}

void kgran_pan(t_kgran *x, double pl, double pm, double ph, double pt) {
	x->panLow = fmax(0, pl);
	x->panMid = fmax(pm, pl);
	x->panHigh = fmin(fmax(ph, pm), 1);
	x->panTight = pt;
}

void kgran_new_grain(t_kgran *x, Grain *grain, double sync){
	//post("New grain!");
	int sr = sys_getsr();
	int head = floor(sync * x->w_len);
	
	double floatShift =  prob(x->grainHeadLow, x->grainHeadMid, x->grainHeadHigh, x->grainHeadTight);
	int idealShift = floor(floatShift * (float)x->w_len);
	//post("2");
	//post("");
	double trans = prob(x->transLow, x->transMid, x->transHigh, x->transTight);
	//post("Passed the trans loop");
	double increment = cpsoct(10.0 + trans) * x->oneover_cpsoct10;
	float offset; // deviation every sample versus the head
	//post("3");
	if (x->w_connected)
		offset = increment - 1; // moving buffer
	else
		offset = increment; // static buffer 
	//post("4");
	float grainDurSamps = (float) prob(x->grainDurLow, x->grainDurMid, x->grainDurHigh, x->grainDurTight) * sr;
	int sampOffset = (int) round(abs(grainDurSamps * offset)); // how many total samples the grain will deviate from the normal buffer movement

	//post("Offset calculated");
	
	if (sampOffset >= x->w_len) // this grain cannot exist with size of the buffer
	{
		error("GRAIN IGNORED, TRANSPOSITION OR DURATION TOO EXTREME");
		return;
	}
	
	int minShift;
	int maxShift;

	if (offset > 0)
	{
		minShift = sampOffset;
		maxShift = x->w_len;
	}
	else
	{
		minShift = 1;
		maxShift = x->w_len - sampOffset;
	}

	if (maxShift == minShift)
	{
		return; // There's a better way to handle this that I'll add at some point...
	}
	
	grain->currTime = head + fmax(fmin(idealShift, maxShift), minShift); // adjust the grain so that it retains its duration within the buffer limitations
	
	//post("Current time set");
	
	
	
	float panR = (float) prob(x->panLow, x->panMid, x->panHigh, x->panTight);
	grain->waveSampInc = increment;
	grain->ampSampInc = ((float)x->w_envlen) / grainDurSamps;
	
	grain->isplaying = true;
	
	grain->ampPhase = 0;
	grain->panR = panR;
	grain->panL = 1 - panR; // separating these in RAM means fewer sample rate calculations
	grain->endTime = grainDurSamps * increment + grain->currTime;
	//post("New grain with start time %f end time %f", grain->currTime, grain->endTime);
}

void kgran_reset_grain_rate(t_kgran *x){
	x->newGrainCounter = (int)round(sys_getsr() * prob(x->grainRateVarLow, x->grainRateVarMid, x->grainRateVarHigh, x->grainRateVarTight));
	//post("Setting grain rate to %d", x->newGrainCounter);
}


// rewrite
void kgran_perform64(t_kgran *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam)
{
	t_double		*in = ins[0];
	t_double		*r_out = outs[0];
	t_double		*l_out = outs[1];
	
	
	int				n = sampleframes;
	float			*b;
	float			*e;
	
	t_buffer_obj	*buffer = buffer_ref_getobject(x->w_buf);
	t_buffer_obj	*env = buffer_ref_getobject(x->w_env);

	b = buffer_locksamples(buffer);
	e = buffer_locksamples(env);
	
	
	double head = 0;
	
	if (!b || !e || !x->running)
	{
		//post("DSP failure");
		goto zero;
	}
	//post("Running DSP loop");
	while (n--){
		if (x->w_connected) // check if the inlet is connected
			head = *in++;
		*r_out = 0; // r_out already stores the input signal for some reason?  so we have to set it to 0
		*l_out = 0;
		for (size_t j = 0; j < MAXGRAINS; j++){
			Grain* currGrain = &x->grains[j];
			if (currGrain->isplaying)
			{
				if ((*currGrain).currTime > currGrain->endTime)
				{
					currGrain->isplaying = false;
					//post("Deleting grain with amp phase %f", currGrain->ampPhase);
				}
				else
				{
					int currTimeInBuffer = (int)floor(currGrain->currTime) % x->w_len;
					float grainAmp = oscili(1, currGrain->ampSampInc, e, x->w_envlen, &((*currGrain).ampPhase));
					float grainOut = grainAmp * b[currTimeInBuffer]; // should include an interpolation option at some point
					currGrain->currTime += currGrain->waveSampInc;
					*l_out += (grainOut * (double)currGrain->panL);
					*r_out += (grainOut * (double)currGrain->panR);
				}
			}
			// this is not an else statement so a grain can be potentially stopped and restarted on the same frame

			if ((x->newGrainCounter <= 0) && !currGrain->isplaying)
			{
				kgran_reset_grain_rate(x);
				if (x->newGrainCounter > 0) // we don't allow two grains to be create on the same frame
					{kgran_new_grain(x, currGrain, head);}
				else
					{x->newGrainCounter = 1;}

			}
		}
		x->newGrainCounter--;
		l_out++;
		r_out++;
	}
	//post("DSP loop finished");
	// if all current grains are occupied, we skip this request for a new grain
	if (x->newGrainCounter <= 0)
	{
		kgran_reset_grain_rate(x);
	}
	
	
	

	buffer_unlocksamples(buffer);
	buffer_unlocksamples(env);
	return;
zero:
	while (n--) {
		*l_out++ = 0.;
		*r_out++ = 0.;
	}
}

// adjust for the appropriate number of inlets and outlets (2 out, one in)
void kgran_dsp64(t_kgran *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags)
{
	x->w_connected = count[1];
	object_method(dsp64, gensym("dsp_add64"), x, kgran_perform64, 0, NULL);
}
