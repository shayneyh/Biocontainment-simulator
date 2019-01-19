/*

This is the model for couple logistic model.

Why we do this?
A safetyguard simulation that can run on the website.

Who is going to see it?
The people who look at our proposal.

What are they expect to see?
In a short time (20s), they can see the difference between two proposals. Visually appeal to people

How would we show this?
two ways:
by coding or by showing a video.

*/

/**********
Parameters
***********/
// Dimmer variable
var dim = 1

// Initial Populaiton
var pop = [190, 10];

// Parameters for the model.
var params = {
	'r1' : 0.2,
	'r2' : 0.1,
	'k1' : 200,
	'k2' : 200,
	'a' : 1,
	'b' : 1};

// The number of steps will be performed.
var steps = 300;

// This is the expected catasrophe step and percentage.
var kill_step = 50;
var kill_percentage = 0.85;

// This is for graphics styles.
var color = ['#DB4105', '#33332D'];

/**********
   Model
***********/

function poisson(expectvalue){
 
    var n = 0, 
	limit = Math.exp(-expectvalue), 
	x = Math.random(); 
 
    while(x > limit){
        n++;
        x *= Math.random();
    }
 
    return n;
};

// Function for couple logisitic model.
function couple_logistic_possison(pop, params, steps) {
	
	/*
	// Nondimensionalize our parameters
	var u1 = pop[0]/params.r1;
	var u2 = pop[1]/params.r2;
	var p = params.r2/params.r1;
	var t = params.r1*cur_t;
	var ab = params.a*params.k2/params.k1;
	var ba = params.b*params.k1/params.k2;

	// Change in population.
	var delta_pop1 = u1*(1 + u1 - ab*u2);
	var delta_pop2 = u2*(1 + u2 - ba*u1);
	*/

	// Change in population.
	var delta_pop1 = params.r1 * pop[0] * (1 - (pop[0] + params.a*pop[1])/params.k1);
	var delta_pop2 = params.r2 * pop[1] * (1 - (pop[1] + params.a*pop[0])/params.k2);
	
	// Call poisson distribution.
	var actual_pop1 = poisson(delta_pop1);
	var actual_pop2 = poisson(delta_pop2);

	// Population growth.
	pop[0] += actual_pop1;
	pop[1] += actual_pop2;

	return [cur_t, pop[0]/dim, pop[1]/dim];

};

// This gives a catastrophe event.
function cata(percentage, pop){
	pop[0] = poisson((1-percentage)*pop[0]);
	pop[1] = poisson((1-percentage)*pop[1]);
}


/**********
 Graphics
***********/

// Draw dots.
function draw(x, y, type, color) {
  	var canvas = document.getElementById(type);
	var c = canvas.getContext('2d');
	c.arc(2*x, 300-y, 0.3, 0, Math.PI*2,true);
	c.strokeStyle = color;
	c.stroke();
}

function show_graphics(result){
	console.log(result);
	document.getElementById("result").innerHTML = "step: " + result[0].toString() +
														"<br>wild type population: " + result[1].toString() + 
														"<br>mutant population: " + result[2].toString();
	draw(result[0], result[1], 'WT', color[0]);
	draw(result[0], result[2], 'MUT', color[1]);
}


// Function that shows the result from logisitic model.
function main(){
	cur_t = 0
	kill_step = poisson(kill_step);
	count_step = 0;
	//do {
	var times = setInterval(function(){
		if (cur_t == steps) clearInterval(times);

		var temp = couple_logistic_possison(pop, params, steps);
		if (count_step == kill_step){
			cata(kill_percentage, pop);
			kill_step = poisson(kill_step);
			count_step = 0;
		}

		show_graphics(temp);

		cur_t += 1;
		count_step += 1;
		return temp;
		}, 35);
	//} while (cur_t < steps)	
};
