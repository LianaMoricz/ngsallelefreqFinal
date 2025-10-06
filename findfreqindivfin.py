import tkinter as tk
from tkinter import filedialog, simpledialog
import os
import json
import webbrowser
import random
import statistics
import matplotlib.pyplot as plt
import numpy as np

def classifyalleles(filename, refsequence):
    cutsitemarker = "GATCGCC"
    cutsiteindex = refsequence.find(cutsitemarker)
    
    crdeletionstart = cutsiteindex + len(cutsitemarker) + 4
    crdeletionend = crdeletionstart + 3
    
    csreads = 0
    crreads = 0
    nhejreads = 0
    otherreads = 0
    topsequences = []
    
    with open(filename, 'r') as f:
        header = next(f) 
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue
                
            alignedseq = parts[0]
            refseqinfile = parts[1]
            reads = int(parts[6])
            percentage = float(parts[7])
            
            topsequences.append({'alignedseq': alignedseq, 'reads': reads, 'percentage': percentage, 'type': 'unknown'})
            #find first character, reference chunk
            refstartinfull = refsequence.find(refseqinfile) 
            if refstartinfull == -1:
                continue
                
            startinaligned = crdeletionstart - refstartinfull
            endinaligned = crdeletionend - refstartinfull
            
            leadinghyphens = 0
            for char in alignedseq:
                if char == '-':
                    leadinghyphens += 1
                else:
                    break

            crdeletionpresent = (startinaligned >= 0 and endinaligned <= len(alignedseq) and alignedseq[startinaligned:endinaligned] == '---')
            
            #cs/nhej/other starting is >5
            totalhyphens = alignedseq.count('-')
            
            #allele type
            if leadinghyphens >= 5:

                remaining = alignedseq[5:]
                if '-' in remaining:  
                    #has other deletions NHEJ
                    topsequences[-1]['type'] = 'nhej'
                    #CS (12bp deletion)
                    #if also CR,other
                    if crdeletionpresent:
                        topsequences[-1]['type'] = 'other'
                        
                        beforeindex = startinaligned - 1
                        afterindex = endinaligned
                        if (beforeindex >= 0 and alignedseq[beforeindex] == '-') or (afterindex < len(alignedseq) and alignedseq[afterindex] == '-'):
                            #nhej if any next/before base is deletion, not actually htr deletion
                            topsequences[-1]['type'] = 'nhej'

                else:
                    #only CS, CS
                    topsequences[-1]['type'] = 'cs'
            elif crdeletionpresent and totalhyphens == 3:
                #CRonly,  CR
                topsequences[-1]['type'] = 'cr'
            else:
                #all other
                topsequences[-1]['type'] = 'other'



            if topsequences[-1]['type'] == 'cs':
                csreads += reads
            if topsequences[-1]['type'] == 'other':
                otherreads += reads
            if topsequences[-1]['type'] == 'nhej':
                nhejreads += reads
            if topsequences[-1]['type'] == 'cr':
                crreads += reads

    
    #top 50
    topsequences.sort(key=lambda x: x['reads'], reverse=True)
    top50 = topsequences[:50]
    
    totalreads = csreads + crreads + otherreads + nhejreads
    totalCORreads = csreads + crreads + nhejreads

    cspercentage = (csreads / totalreads) * 100 if totalreads > 0 else 0
    crpercentage = (crreads / totalreads) * 100 if totalreads > 0 else 0
    otherpercentage = (otherreads / totalreads) * 100 if totalreads else 0 
    nhejpercentage = (nhejreads / totalreads) * 100 if totalreads else 0 

    cspercentagecorr = (csreads /totalCORreads) * 100 if totalCORreads > 0 else 0
    crpercentagecorr = (crreads / totalCORreads) * 100 if totalCORreads > 0 else 0
    nhejpercentagecorr = (nhejreads / totalCORreads) * 100 if totalCORreads else 0 


    HTRpercentage = (crreads - (csreads + nhejreads)) / ( crreads + csreads + nhejreads) * 100
    # if HTRpercentage < 0:
    #     HTRpercentage = 0
    
    return {
        'filename': os.path.basename(filename),
        'csreads': csreads, #good
        'cspercentage': cspercentage, #good
        'crreads': crreads, #good
        'crpercentage': crpercentage, #good
        'otherreads': otherreads, #good
        'otherpercentage': otherpercentage, #good
        'totalreads': totalreads, #good
        'topsequences': top50, #good
        'nhejreads': nhejreads, #good
        'nhejpercentage': nhejpercentage, #good
        'HTRpercentage': HTRpercentage, #good
        'cspercentagecorr': cspercentagecorr, #good
        'crpercentagecorr': crpercentagecorr, #good
        'nhejpercentagecorr': nhejpercentagecorr, #good
        'totalCORreads': totalCORreads #good
        
    }

def generatehtmlreport(genotyperesults, outputfile):
    allresults = []
    for group in genotyperesults:
        for result in group['samples']:
            resultWithGroup = result.copy()
            resultWithGroup['group'] = group['name']
            allresults.append(resultWithGroup)

    #bar chart with error bars
    groupNames = [group['name'] for group in genotyperesults]
    
    #means and stdev for each group/type
    barChartData = {}
    metrics = ['cspercentagecorr', 'nhejpercentagecorr', 'crpercentagecorr', 'HTRpercentage']
    metricLabels = ['CS Corrected %', 'NHEJ Corrected %', 'CR Corrected %', 'HTR %']
    
    for group in genotyperesults:
        groupName = group['name']
        barChartData[groupName] = {}
        
        for metric in metrics:
            values = [s[metric] for s in group['samples']]
            if values:
                barChartData[groupName][metric] = {
                    'mean': statistics.mean(values),
                    'stdev': statistics.stdev(values) if len(values) > 1 else 0
                }
            else:
                barChartData[groupName][metric] = {
                    'mean': 0,
                    'stdev': 0
                }

    htmlcontent = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Allele Frequency Report from CRISPResso Individual </title>
        <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
        <style>
            body {{
                font-family: Lucida Console, monospace;
                margin: 0.05;
                padding: 16px;
                background-color: #f5f5f5;
            }}
            .container {{
                max-width: 1200px;
                margin: 0 auto;
                background-color: white;
                padding: 20px;
                border-radius: 8px;
                box-shadow: 0 2px 10px rgba(0, 0, 0, 0.1);
            }}
            h1, h2, h3 {{
                color: #333;
            }}
            .chart-container {{
                margin: 30px 0;
                position: relative;
                height: 500px;
            }}
            table {{
                width: 100%;
                border-collapse: collapse;
                margin: 20px 0;
            }}
            th, td {{
                padding: 12px;
                text-align: left;
                border-bottom: 1px solid #ddd;
            }}
            th {{
                background-color: #f2f2f2;
            }}
            tr:hover {{
                background-color: #f5f5f5;
            }}
            .summary {{
                background-color: #d3d3d3;
                padding: 15px;
                border-radius: 5px;
                margin: 20px 0;
            }}
            .sequence-box {{
                font-family: monospace;
                background-color: #f9f9f9;
                padding: 10px;
                border-radius: 5px;
                margin: 10px 0;
                overflow-x: auto;
            }}
            .cs-seq {{
                color: #415b76;
            }}
            .cr-seq {{
                color: #c0392b;
            }}
            .nhej-seq {{
                color: #27ae60;
            }}
            .htr-seq {{
                color: #f39c12;
            }}
            .other-seq {{
                color: #7f8c8d;
            }}
            .sample-section {{
                margin: 30px 0;
                padding: 20px;
                border: 1px solid #ddd;
                border-radius: 5px;
            }}
            .collapsible {{
                background-color: #eee;
                color: #444;
                cursor: pointer;
                padding: 18px;
                width: 100%;
                border: none;
                text-align: left;
                outline: none;
                font-size: 15px;
                margin: 10px 0;
            }}
            .active, .collapsible:hover {{
                background-color: #ccc;
            }}
            .content {{
                padding: 0 18px;
                display: none;
                overflow: hidden;
                background-color: #f9f9f9;
            }}
        </style>
    </head>
    <body>
        <div class="container">
            <h1>Allele Frequency Report from CRISPResso Individual</h1>
            
            <div class="Check:">
                <h2></h2>
                <p>{len(allresults)} files across {len(genotyperesults)} groups.</p>
            </div>

            <h2>Corrected Percentages by Genotype</h2>
            <div class="chart-container">
                <canvas id="barChart"></canvas>
            </div>

            <h2>Percentages</h2>
            <table>
                <thead>
                    <tr>
                        <th>Group</th>
                        <th>Filename</th>
                        <th>Total Reads</th>
                        <th>Total Corrected Reads</th>
                        <th>CS Reads</th>
                        <th>CS %</th>
                        <th>CS Corrected %</th>
                        <th>NHEJ Reads</th>
                        <th>NHEJ %</th>
                        <th>NHEJ Corrected %</th>
                        <th>CR Reads</th>
                        <th>CR %</th>
                        <th>CR Corrected %</th>
                        <th>HTR %</th>
                        <th>Other Reads</th>
                        <th>Other %</th>
                    </tr>
                </thead>
                <tbody>
    """   

    for result in allresults:
        htmlcontent += f"""
                    <tr>
                        <td>{result['group']}</td>
                        <td>{result['filename']}</td>
                        <td>{result['totalreads']}</td>
                        <td>{result['totalCORreads']}</td>
                        <td>{result['csreads']}</td>
                        <td>{result['cspercentage']:.2f}%</td>
                        <td>{result['cspercentagecorr']:.2f}%</td>
                        <td>{result['nhejreads']}</td>
                        <td>{result['nhejpercentage']:.2f}%</td>
                        <td>{result['nhejpercentagecorr']:.2f}%</td>
                        <td>{result['crreads']}</td>
                        <td>{result['crpercentage']:.2f}%</td>
                        <td>{result['crpercentagecorr']:.2f}%</td>
                        <td>{result['HTRpercentage']:.2f}%</td>
                        <td>{result['otherreads']}</td>
                        <td>{result['otherpercentage']:.2f}%</td>
                    </tr>
        """

    htmlcontent += """
                </tbody>
            </table>
            
            <h2>Top 50 Sequences by Sample</h2>
    """
        
    # Top sequences for each sample
    for result in allresults:
        htmlcontent += f"""
            <button class="collapsible">{result['group']} - {result['filename']} - Top 50 Sequences</button>
            <div class="content">
                <table>
                    <thead>
                        <tr>
                            <th>Type</th>
                            <th>Sequence</th>
                            <th>Reads</th>
                            <th>Percentage</th>
                        </tr>
                    </thead>
                    <tbody>
        """

        for seq in result['topsequences']:
            seqclass = f"{seq['type']}-seq"
            htmlcontent += f"""
                        <tr>
                            <td><span class="{seqclass}">{seq['type'].upper()}</span></td>
                            <td class="sequence-box {seqclass}">{seq['alignedseq']}</td>
                            <td>{seq['reads']}</td>
                            <td>{seq['percentage']:.2f}%</td>
                        </tr>
            """

        htmlcontent += """
                    </tbody>
                </table>
            </div>
        """  

    htmlcontent += f"""
            <script>
                const barChartData = {json.dumps(barChartData)};
                const groupNames = {json.dumps(groupNames)};
                const metrics = {json.dumps(metrics)};
                const metricLabels = {json.dumps(metricLabels)};
                
                // Define colors for each metric type
                const metricColors = [
                    'rgba(65, 91, 118, 0.7)',    // CS - blue
                    'rgba(39, 174, 96, 0.7)',    // NHEJ - green
                    'rgba(192, 57, 43, 0.7)',    // CR - red
                    'rgba(243, 156, 18, 0.7)'    // HTR - orange
                ];
                
                const borderColors = [
                    'rgba(65, 91, 118, 1)',
                    'rgba(39, 174, 96, 1)',
                    'rgba(192, 57, 43, 1)',
                    'rgba(243, 156, 18, 1)'
                ];
                
                // Prepare datasets for bar chart
                const datasets = [];
                
                metrics.forEach((metric, metricIndex) => {{
                    const data = [];
                    const errorBars = [];
                    
                    groupNames.forEach(groupName => {{
                        if (barChartData[groupName] && barChartData[groupName][metric]) {{
                            data.push(barChartData[groupName][metric].mean);
                            errorBars.push(barChartData[groupName][metric].stdev);
                        }} else {{
                            data.push(0);
                            errorBars.push(0);
                        }}
                    }});
                    
                    datasets.push({{
                        label: metricLabels[metricIndex],
                        data: data,
                        backgroundColor: metricColors[metricIndex],
                        borderColor: borderColors[metricIndex],
                        borderWidth: 1,
                        errorBars: errorBars
                    }});
                }});
                
                // Create bar chart
                const barCtx = document.getElementById('barChart').getContext('2d');
                new Chart(barCtx, {{
                    type: 'bar',
                    data: {{
                        labels: groupNames,
                        datasets: datasets
                    }},
                    options: {{
                        responsive: true,
                        maintainAspectRatio: false,
                        scales: {{
                            x: {{
                                title: {{
                                    display: true,
                                    text: 'Genotype Groups'
                                }}
                            }},
                            y: {{
                                beginAtZero: true,
                                max: 100,
                                title: {{
                                    display: true,
                                    text: 'Percentage'
                                }}
                            }}
                        }},
                        plugins: {{
                            legend: {{
                                position: 'top'
                            }},
                            tooltip: {{
                                callbacks: {{
                                    label: function(context) {{
                                        const datasetIndex = context.datasetIndex;
                                        const groupIndex = context.dataIndex;
                                        const mean = context.parsed.y;
                                        const stdev = datasets[datasetIndex].errorBars[groupIndex];
                                        return metricLabels[datasetIndex] + ': ' + mean.toFixed(2) + '% ± ' + stdev.toFixed(2) + '%';
                                    }}
                                }}
                            }}
                        }}
                    }}
                }});
                
                // Collapsible sections
                var coll = document.getElementsByClassName("collapsible");
                var i;
                
                for (i = 0; i < coll.length; i++) {{
                    coll[i].addEventListener("click", function() {{
                        this.classList.toggle("active");
                        var content = this.nextElementSibling;
                        if (content.style.display === "block") {{
                            content.style.display = "none";
                        }} else {{
                            content.style.display = "block";
                        }}
                    }});
                }}
            </script>
        </body>
    </html>
    """

    with open(outputfile, 'w') as f:
        f.write(htmlcontent)
    return outputfile

def generatemplchart(genotyperesults):
    
    groupnames = [group['name'] for group in genotyperesults]
    metrics = ['cspercentagecorr', 'nhejpercentagecorr', 'crpercentagecorr', 'HTRpercentage']
    metriclabels = ['CS Corrected %', 'NHEJ Corrected %', 'CR Corrected %', 'HTR %']
    colors = ['#415b76', '#27ae60', '#c0392b', '#f39c12']
    
    means = {metric: [] for metric in metrics}
    stds = {metric: [] for metric in metrics}
    
    for group in genotyperesults:
        for metric in metrics:
            values = [sample[metric] for sample in group['samples']]
            if values:
                means[metric].append(statistics.mean(values))
                stds[metric].append(statistics.stdev(values) if len(values) > 1 else 0)
            else:
                means[metric].append(0)
                stds[metric].append(0)
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    ngroups = len(groupnames)
    nmetrics = len(metrics)
    barwidth = 0.2
    index = np.arange(ngroups)
    
    bars = []
    for i, metric in enumerate(metrics):
        barpositions = index + i * barwidth
        bar = ax.bar(barpositions, means[metric], barwidth, 
                    yerr=stds[metric], capsize=5, label=metriclabels[i],
                    color=colors[i], alpha=0.7, edgecolor='black', linewidth=0.5)
        bars.append(bar)
        
        for j, (mean_val, std_val) in enumerate(zip(means[metric], stds[metric])):
            ax.text(barpositions[j], mean_val + std_val + 1, 
                   f'{mean_val:.1f}±{std_val:.1f}', 
                   ha='center', va='bottom', fontsize=9)
    
    ax.set_xlabel('Genotype Groups', fontsize=12, fontweight='bold')
    ax.set_ylabel('Percentage', fontsize=12, fontweight='bold')
    ax.set_title('Corrected Percentages by Genotype', fontsize=14, fontweight='bold')
    ax.set_xticks(index + barwidth * (nmetrics - 1) / 2)
    ax.set_xticklabels(groupnames)
    ax.set_ylim(0, 100)
    ax.legend()
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.show()


def main():
    root = tk.Tk()
    root.withdraw()
    #window will open
    #refsequence defined for for finding positions afterwards.
    refsequence = "TTGCGGCGTGGCCTATCCGGGCGAACTTTTGGCCGTGATGGGCAGTTCCGGTGCCGGAAAGACGACCCTGCTGAATGCCCTTGCCTTTCGATCGCCGCAGGGCATCCAAGTATCGCCATCCGGGATGCGACTGCTCAATGGCCAACCTGTGGACGCCAAGGAGATGCAGGCCAGGTGCGCCTATGTCCAGCAGGATGACCTCTTTATCGGCTCCCTAACGGCCAGGGAACACCTGATTTTCCAAGCCATGGTGCGGATGCCACGACATCTGACCTATCGGCAGCGAGTGGCCCGCGTGGATCAGGTGATCCAGGAGCTTTCGCTCAGCAAATGTCAGCACACGATCATCGGTGTGCCCGGCAGGGTGAAAGGTCTGTCCGGCGGAGAAAGG"
    
    genotyperesults = []

    #while window is open:
    while True:
        genotypename = simpledialog.askstring("Group", "Name/Cancel:")
        if genotypename is None:
            break
        filenames = filedialog.askopenfilenames(title=f"select .txt Files for {genotypename}", filetypes=[("Text files", "*.txt")])
        if not filenames:
            continue
        
        genotypesamplelist = []
        
        for filename in filenames:
            print(f"file: {os.path.basename(filename)} for GROUP = {genotypename}")
            eachresult = classifyalleles(filename, refsequence)
            genotypesamplelist.append(eachresult)
            print(f" CS reads: {eachresult['csreads']}, CS %: {eachresult['cspercentage']:.2f}%")
            print(f" CR reads: {eachresult['crreads']}, CR %: {eachresult['crpercentage']:.2f}%")
            print(f" NHEJ reads: {eachresult['nhejreads']}, NHEJ %: {eachresult['nhejpercentage']:.2f}%")
            print(f" Other reads: {eachresult['otherreads']}, Other %: {eachresult['otherpercentage']:.2f}%")
            print(f" Total reads: {eachresult['totalreads']}")
            print(f" HTR percentage: {eachresult['HTRpercentage']:.2f}%")
            print(f" Corr CS %: {eachresult['cspercentagecorr']:.2f}%")
            print(f" Corr CR %: {eachresult['crpercentagecorr']:.2f}%")
            print(f" Corr NHEJ %: {eachresult['nhejpercentagecorr']:.2f}%")
            print(f" Total Corr Reads: {eachresult['totalCORreads']:.2f}")
        
        genotyperesults.append({
            'name': genotypename,
            'samples': genotypesamplelist
        })
    
    if not genotyperesults:
        return
    
    generatemplchart(genotyperesults)


    htmloutput = "htmloutput.html"
    reportpath = generatehtmlreport(genotyperesults, htmloutput)
    print(f"\nPath =  {reportpath}")
    webbrowser.open('file://' + os.path.realpath(reportpath))

if __name__ == "__main__":
    main()