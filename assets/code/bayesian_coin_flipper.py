"""
The Bayesian Coin-Flipper
--------------------------------------------------------------------------------

Author: Griffin Chure 
Last Modified: September 23, 2019
License: MIT

Description
--------------------------------------------------------------------------------
This Python script generates an interactive tool that illustrates the principle
of using Bayes' theorem to estimate the bias of a coin. 
"""
import numpy as np
import bokeh.io
import bokeh.plotting
from bokeh.models import *
import bokeh.layouts
colors = {'black':'#444147', 'purple': '#7E59A2', 'orange':'#E39943', 
                'light_purple':'#A17DB8', 'light_orange':'#EEBA7F'}
bokeh.io.output_file('bayesian_coin_flipper.html')


# Define sliders. 
initialize = Button(label='Click to generate new coinflips')
coin_flips = Slider(title='Display first log\u2081\u2080 coin flips', 
                    start=0, end=6, step=0.001, value=1, 
                    bar_color=colors['light_purple'])
coin_bias = Slider(title='Coin bias', start=0.001, end=1, step=0.001, value=0.5,
                    bar_color=colors['light_purple'])
prior_mu = Slider(title='Prior µ', start=0, end=1, step=0.001, value=0.5,
                    bar_color=colors['light_orange'])
prior_sig = Slider(title='Prior σ', start=0.0001, end=5, step=0.0001, value=0.1,
                    bar_color=colors['light_orange'])

# Set up a plot showing the bias and calculated bias as the number of flips
post_ax = bokeh.plotting.figure(width=550, height=300, 
                            y_axis_label='\u221d probability', 
                            x_axis_label='coin bias',
                            x_range=[0, 1], y_range=[0, 1.05])

# Format the y axis ticks. 
post_ax.ygrid.visible = False
post_ax.yaxis.ticker = [-10, 10]
post_ax.title.text_font_style = 'italic'
post_ax.title.background_fill_alpha = 0
post_ax.title.align = 'left'
post_ax.title.text = "Click button to initialize flips"


# Define a data source for showing the posterior and prior
n_points = 500
post_source = ColumnDataSource({'probability': np.linspace(0.001, 0.999, n_points), 
                                'posterior': np.zeros(n_points), 
                                'prior':np.zeros(n_points)})
# Define a source to display the true coin bias
bias_source = ColumnDataSource({'bias':[0.5]})
flip_source = ColumnDataSource({'flips':[]})


# Define the js for changing the coin bias
bias_js = CustomJS(args={'bias_source':bias_source, 'coin_bias':coin_bias},
                    code="""
                    bias_source.data['bias'] = [coin_bias.value]; 
                    bias_source.change.emit();""") 

# Define JS code to view a subset of the coin flips and compute the number of 
# heads
viewer_cb = """
    // Filters the output of flipper based on the number of flips to observe. 
    var n_displayed = Math.floor(Math.pow(10, flip_viewer.value));

    // compute the number of heads in the newly sliced array
    var displayed_flips = flip_source.data['flips'].slice(0, n_displayed);
    var n_flips = n_displayed;
    var n_heads = displayed_flips.reduce((v1, v2) => v1 + v2);
    
    // Update the title to reflect the number of displayed coin flips.
    title.text = "Displaying " + n_displayed + " coin flips"
"""

# Define JS code to generate one million coin flips with the provided bias.
flipper_cb = """
// Generate the coinflips based on the entered values   
var bias = coin_bias.value;

// *Always* compute 1E6 flips
var flips = new Array(1000000).fill().map(() => Math.random() < bias);
flip_source.data['flips'] = flips;
flip_source.change.emit();
"""

# Define the JS code to perform the inference. Note that factorials are 
# approximated by using a gamma function.
inference_cb = """
// Define a gamma function using the Lanczos approximation
// See https://github.com/substack/gamma.js for node module
 var g_ln = 607/128;
 var p_ln = [0.99999999999999709182,57.156235665862923517,-59.597960355475491248,
            14.136097974741747174,-0.49191381609762019978,0.33994649984811888699e-4,
            0.46523628927048575665e-4,-0.98374475304879564677e-4,0.15808870322491248884e-3,
            -0.21026444172410488319e-3,0.21743961811521264320e-3,-0.16431810653676389022e-3,
            0.84418223983852743293e-4, -0.26190838401581408670e-4,0.36899182659531622704e-5];
 // Spouge approximation (suitable for large arguments)
function lngamma(z) {  
      if (z < 0) {return Number('0/0')};
      var x = p_ln[0];
      for(var i = p_ln.length - 1; i > 0; --i) { 
            x += p_ln[i] / (z + i)
        };
      var t = z + g_ln + 0.5;
      return .5*Math.log(2*Math.PI)+(z+.5)*Math.log(t)-t+Math.log(x)-Math.log(z);
  }


// Define the log likelihood
function logLike(n, N, p) {
    var binomCoeff = lngamma(N + 1) - lngamma(n + 1) - lngamma(N - n + 1)
    var prob = n * Math.log(p) +  (N - n) * Math.log(1 - p)
    return binomCoeff + prob;
}

// Define the log prior
function logPriorNormal(x) { 
    var prefactor = -0.5 * Math.log(Math.sqrt(2 * Math.PI * Math.pow(prior_sig.value, 2)));
    var expon = - Math.pow(x - prior_mu.value, 2) / (Math.pow(prior_sig.value, 2));
    return prefactor + expon
}

function logSumExp(vals) {
    var maxVal = Math.max.apply(null, vals);
    var sum_value = 0
    for (var i = 0; i < vals.length; i++) {
        sum_value += Math.exp(vals[i] - maxVal);
    } 
   return maxVal + Math.log(sum_value);
}

// Evaluate the posterior
var log_posterior = [];
var prior = []
for (var i = 0; i < post_source.data['probability'].length; i++) {
    var prob = post_source.data['probability'][i]
    var log_likelihood = logLike(n_heads, n_flips, prob);

    log_prior = logPriorNormal(prob);
    log_posterior[i] = log_likelihood + log_prior;
    prior[i] = Math.exp(log_prior);
}

// Use the logsumexp trick to normalize the posterior
var summed = logSumExp(log_posterior);
var posterior= []
for (var i = 0 ; i < log_posterior.length; i++) { 
   posterior[i] = Math.exp(log_posterior[i] - summed);
}

// Rescale everything to be on the same meaningless scale ([0, 1]).
var maxPrior = Math.max(...prior)
var maxPosterior = Math.max(...posterior)
var priorNorm = []
var posteriorNorm = []
for (var i = 0; i < post_source.data['probability'].length; i++) {
    priorNorm[i] = prior[i] / maxPrior;
    posteriorNorm[i] = posterior[i] / maxPosterior;
}

// Update the source data.
post_source.data['posterior'] = posteriorNorm;
post_source.data['prior'] = priorNorm;
post_source.change.emit();
"""

# Define the callbacks
args = {'coin_bias':coin_bias, 
        'prior_mu':prior_mu, 
        'prior_sig':prior_sig, 
        'flip_viewer':coin_flips,
        'post_source':post_source,
        'flip_source':flip_source,
        'title': post_ax.title}
flipper = CustomJS(args=args, code=flipper_cb + viewer_cb + inference_cb) 
viewer = CustomJS(args=args, code=viewer_cb + inference_cb)

# Assign the callbacks to the inputs
initialize.js_on_click(flipper)
coin_bias.js_on_change('value', bias_js)
coin_flips.js_on_change('value', viewer)
prior_sig.js_on_change('value', viewer)
prior_mu.js_on_change('value', viewer)


# Define the plot layout
box1 = WidgetBox(coin_flips, coin_bias)
box2 = WidgetBox(prior_mu, prior_sig)
row1 = bokeh.layouts.row(box1, box2)
lay = bokeh.layouts.column(initialize, row1, post_ax)


# Populate the axes
post_ax.line(x='probability', y='posterior', source=post_source,
            color=colors['purple'], line_width=2, legend='posterior probability')
post_ax.line(x='probability', y='prior', source=post_source, 
            color=colors['orange'], line_width=2, legend='prior probability')
post_ax.ray(x='bias', y=0, source=bias_source, 
            length=0, angle=np.pi/2, line_width=2, color=colors['black'],
            legend='true coin bias')

# Set the theme details
theme_json =  {
    'attrs' : {
        'Figure' : {
            'background_fill_color': '#EEEEEE',
        },
        'Axis': {
            'axis_line_color': 'slategray',
            'major_tick_line_color': None,
            'minor_tick_line_color': None,
        },
        'Legend': {
            'border_line_color': 'slategray',
            'background_fill_color': '#EEEEEE',
            'border_line_width': 0.75,
            'background_fill_alpha': 0.75,
        },
        'Grid': {
            'grid_line_color': '#FFFFFF',
            'grid_line_width': 0.75,
        },
        'Text': {
            'text_font_style': 'italic',
            'text_font': 'Arial', 
            'text_font_size':10,
        },
        'Title': {
            'text_color': '#3c3c3c',
            'text_font': 'Arial',
            'offset': 2,
         }
    }
}

        
theme = bokeh.themes.Theme(json=theme_json)
bokeh.io.curdoc().theme = theme
bokeh.io.save(lay)

