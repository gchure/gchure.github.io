"""
Gillespie Simulation of a Constitutive Promoter
--------------------------------------------------------------------------------
Author: Griffin Chure
Last Modified: September 24, 2019
License: MIT

Description
--------------------------------------------------------------------------------
This script generates a stand-alone interactive tool that performs a Gillespie
simulation of a constitutively expressing promoter, keeping track of the number
of mRNAs produced as a function of time. It is not meant to be scientifically
useful, but is more for illustrative purposes.

Notes
--------------------------------------------------------------------------------
This script has a bug somewhere where the *final* mRNA number 
(after n_traj events) is *always* even. I'm not sure where this bug is lurking.
Again, this tool is for illustrative purposes and should therefore not be used 
as the basis for any *real* Gillespie simulation (not that you would do it in 
JavaScript anyway!)
"""

import numpy as np
import bokeh.io
import bokeh.plotting
from bokeh.models import *
bokeh.io.output_file('gillespie.html')
colors = {'black':'#444147', 'purple': '#7E59A2', 'orange':'#E39943', 
                'light_purple':'#A17DB8', 'light_orange':'#EEBA7F'}
                
# Generate some data sources
source = ColumnDataSource({'mRNAs':[[]], 'times':[[]]})
theo_source = ColumnDataSource({'mRNA':[],  'time':[]})

# Define interactivity inputs
n_traj = Slider(title='number of trajectories', start=10, end=1000, step=5, value=500)
n_sims = Slider(title='number of simulations', start=1, end=1000, step=10, value=100)
r = Slider(title='mRNA production rate (r) [per min.]', 
            start=1, end=50, step=1, value=10)
m_0 = Slider(title='inital number of mRNAs', start=0, end=100, step=1, value=0)
gamma = Slider(title='mRNA degradation rate (γ) [per min.]',
            start=1/30, end=3, step=0.1, value=1/3)

# Define the "run" button
run = Button(label='click to run Gillespie simulation')

# load the gillespie functions
with open('gillespie.js', 'r') as file:
    gillespie_fn = file.read()

# Define the actual gillespie callback
cb_code = """
function randomExponential(rate) {
    rate = rate || 1; 
    return -Math.log(Math.random()) / rate;
}

function runGillespie(num_traj, prod_rate, deg_rate, m_0) { 
    var m_t = [m_0];
    var time = [0];

    // Loop through each trajectory and run the simulation.
    for (var i = 1; i <= num_traj; i++) { 
        // Calculate the propensity for each reaction
        var lam_deg = m_t[i - 1] * deg_rate;

        // Sum the propensities for each reaction
        var sum_propen = prod_rate + lam_deg;

        // Keep track of the current time. 
        time.push(time[i - 1] + randomExponential(sum_propen));

        // Draw a random number and decide which reaction should occur. 
        var draw = Math.random();
        if (draw < (prod_rate / sum_propen)) { 
            // Produce an mRNA
            m_t.push(m_t[i - 1] + 1);}
        else {
            // Degrade an mRNA
            m_t.push(m_t[i - 1] - 1);}
        }
    return { mRNA : m_t, time: time}
  }
// Instantiate empty arrays for the actual gillespie
var mRNAs = [];
var time = [];

// Keep track of what the maximum  and minimum mRNA counts are at the end.

// Iterate through each simulation.
for (var i = 0; i < n_sims.value; i++) { 
    // Run the Gillespie
    results = runGillespie(n_traj.value, r.value, gamma.value, m_0.value);
    mRNAs.push(results['mRNA']);
    time.push(results['time']);
}


// Compute the theoretical result given the inputs using Euler integration
var theo_times = [0];
var theo_mRNAs = [m_0.value];
var dt = 1 / 30; // Hard-code a time step smaller than the probabilities 

for (var i = 1; i < 2 * n_traj.value; i++) {
    theo_mRNAs[i] = theo_mRNAs[i - 1] + r.value * dt - gamma.value * theo_mRNAs[i-1] * dt;
    theo_times[i] = theo_times[i - 1] + dt;
} 

// Update the theory source. 
theo_source.data['mRNA'] = theo_mRNAs;
theo_source.data['time'] = theo_times;

// Update the gillespie source.
source.data['mRNAs'] = mRNAs;
source.data['times'] = time;

// Update the distribution source
source.change.emit();
theo_source.change.emit()
"""

# Define the arguments
args = {'n_traj':n_traj, 'r':r, 'gamma':gamma, 
        'm_0':m_0, 'n_sims':n_sims, 'source':source,
        'theo_source':theo_source}
cb = CustomJS(args=args, code = cb_code);

# Assign the callback to the button. 
run.js_on_click(cb)

# Instantiate the figure canvas
ax = bokeh.plotting.figure(width=400, height=300, 
                        x_axis_label='time [min.]', 
                        y_axis_label='number of mRNAs')


ax.multi_line('times', 'mRNAs', source=source, color=colors['purple'], 
             line_width=0.5, alpha=0.5, level='underlay', 
             legend='simulation')
ax.line('time', 'mRNA', source=theo_source, color=colors['orange'],
         legend='theory', line_width=2)

# Style the legend
ax.legend.orientation = 'horizontal'
# Define the layout            
col1 = bokeh.layouts.column(run, n_traj, n_sims, r, gamma, m_0)
lay = bokeh.layouts.row(col1, ax)


# Define the theme and save
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