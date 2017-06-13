
function drawBarChartH(dataset, w, h, hgap, barPadding) {
  var barWidth = w / dataset.length;
  var fmt = d3.format(".1f");
  var runtimeScale = d3.scale.linear().domain([0, 115]).range([0, h - hgap]);

  var svg = d3.select("body")
    .append("svg")
    .attr("width", w)
    .attr("height", h);

  function yStart(runtime) { return h - hgap - runtimeScale(runtime); }

  var comb_insert_rects = svg.append("g").selectAll("rect").data(dataset).enter().append("rect").attr({
    x: function (d, i) { return i * barWidth; },
    y: function (d) { return yStart(d.comb_insert_time); },
    width: (barWidth - barPadding) / 2,
    height: function (d, i) { return runtimeScale(d.comb_insert_time); },
    fill: function(d) { return "blue"; },
  });
  var comb_query_rects = svg.append("g").selectAll("rect").data(dataset).enter().append("rect").attr({
    x: function (d, i) { return i * barWidth; },
    y: function (d, i) { return yStart(d.comb_insert_time + d.comb_query_time); },
    width: (barWidth - barPadding) / 2,
    height: function (d) { return runtimeScale(d.comb_query_time); },
    fill: function(d) { return "cyan"; },
  });
  var comb_texts = svg.append("g").selectAll("text").data(dataset).enter().append("text")
    .text(function (d) { return fmt(d.comb_insert_time + d.comb_query_time); })
    .attr({
      x: function(d, i) { return i * barWidth + (barWidth - barPadding) / 4; },
      y: function(d) { return yStart(d.comb_insert_time + d.comb_query_time) - 5; },
      'text-anchor': 'middle',
      color: 'blue',
      fill: 'blue',
      'font-family': 'Arial',
      'font-size': '14px',
      'font-weight': 'bold',
    });

  var sort_time_rects = svg.append("g").selectAll("rect").data(dataset).enter().append("rect").attr({
    x: function (d, i) { return i * barWidth + (barWidth - barPadding) / 2; },
    y: function (d, i) { return yStart(d.sort_time); },
    width: (barWidth - barPadding) / 2,
    height: function (d, i) { return runtimeScale(d.sort_time); },
    fill: function(d) { return "red"; },
  });
  var sort_query_time_rects = svg.append("g").selectAll("rect").data(dataset).enter().append("rect").attr({
    x: function (d, i) { return i * barWidth + (barWidth - barPadding) / 2; },
    y: function (d, i) { return yStart(d.sort_time + d.sort_query_time); },
    width: (barWidth - barPadding) / 2,
    height: function (d, i) { return runtimeScale(d.sort_query_time); },
    fill: function(d) { return "magenta"; },
  });
  var sort_texts = svg.append("g").selectAll("text").data(dataset).enter().append("text")
    .text(function (d) { return fmt(d.sort_time + d.sort_query_time); })
    .attr({
      x: function(d, i) { return i * barWidth + (barWidth - barPadding) * 3 / 4; },
      y: function(d) { return yStart(d.sort_time + d.sort_query_time) - 5; },
      'text-anchor': 'middle',
      fill: 'red',
      color: 'red',
      'font-family': 'Arial',
      'font-size': '14px',
      'font-weight': 'bold',
    });

  var query_texts = svg.append("g").selectAll("text").data(dataset).enter().append("text")
    .text(function (d) {
      var p = Math.floor(Math.log(d.q) / Math.log(10));
      if (p==0) return 'Q = 1';
      if (p==1) return 'Q = 10';
      return 'Q = 10^' + p;
    })
    .attr({
      x: function(d, i) { return i * barWidth + (barWidth - barPadding) / 2; },
      y: function(d) { return h - hgap + 15; },
      'text-anchor': 'middle',
      'font-family': 'Arial',
      'font-size': '14px',
      // 'font-weight': 'bold',
    });

  svg.append("g").append("text")
    .text("N = 10^8 (initially unordered), Q = Number of random point queries")
    .attr({
      x: function(d) { return w / 2; },
      y: function(d) { return h - 10; },
      'text-anchor': 'middle',
      'font-family': 'Arial',
      'font-size': '14px',
      'font-weight': 'bold',
    });
}


function drawBarChartV(dataset, w, h, hgap, barPadding) {
  // var barWidth = h / dataset.length;
  // var fmt = d3.format(".1f");
  // var runtimeScale = d3.scale.linear().domain([0, 115]).range([0, w - hgap]);
  // var svg = d3.select("body").append("svg").attr("width", w).attr("height", h);

  // function xStart(runtime) { return hgap + runtimeScale(runtime); }

  // var comb_insert_rects = svg.append("g").selectAll("rect").data(dataset).enter().append("rect").attr({
  //   x: function (d) { return xStart(0); },
  //   y: function (d, i) { return i * barWidth; },
  //   width: function (d, i) { return runtimeScale(d.comb_insert_time); },
  //   height: (barWidth - barPadding) / 2,
  //   fill: function(d) { return "blue"; },
  // });
  // var comb_query_rects = svg.append("g").selectAll("rect").data(dataset).enter().append("rect").attr({
  //   x: function (d) { return xStart(d.comb_insert_time); },
  //   y: function (d, i) { return i * barWidth; },
  //   width: function (d) { return runtimeScale(d.comb_query_time); },
  //   height: (barWidth - barPadding) / 2,
  //   fill: function(d) { return "cyan"; },
  // });
  // var comb_texts = svg.append("g").selectAll("text").data(dataset).enter().append("text")
  //   .text(function (d) { return fmt(d.comb_insert_time + d.comb_query_time); })
  //   .attr({
  //     x: function(d) { return xStart(d.comb_insert_time + d.comb_query_time) + 5; },
  //     y: function(d, i) { return i * barWidth + (barWidth - barPadding) / 4; },
  //     color: 'blue',
  //     fill: 'blue',
  //     'font-family': 'Arial',
  //     'font-size': '14px',
  //     'font-weight': 'bold',
  //   });

  // var sort_time_rects = svg.append("g").selectAll("rect").data(dataset).enter().append("rect").attr({
  //   x: function (d, i) { return i * barWidth + (barWidth - barPadding) / 2; },
  //   y: function (d, i) { return xStart(d.sort_time); },
  //   width: (barWidth - barPadding) / 2,
  //   height: function (d, i) { return runtimeScale(d.sort_time); },
  //   fill: function(d) { return "red"; },
  // });
  // var sort_query_time_rects = svg.append("g").selectAll("rect").data(dataset).enter().append("rect").attr({
  //   x: function (d, i) { return i * barWidth + (barWidth - barPadding) / 2; },
  //   y: function (d, i) { return xStart(d.sort_time + d.sort_query_time); },
  //   width: (barWidth - barPadding) / 2,
  //   height: function (d, i) { return runtimeScale(d.sort_query_time); },
  //   fill: function(d) { return "magenta"; },
  // });
  // var sort_texts = svg.append("g").selectAll("text").data(dataset).enter().append("text")
  //   .text(function (d) { return fmt(d.sort_time + d.sort_query_time); })
  //   .attr({
  //     x: function(d, i) { return i * barWidth + (barWidth - barPadding) * 3 / 4; },
  //     y: function(d) { return xStart(d.sort_time + d.sort_query_time) - 5; },
  //     'text-anchor': 'middle',
  //     fill: 'red',
  //     color: 'red',
  //     'font-family': 'Arial',
  //     'font-size': '14px',
  //     'font-weight': 'bold',
  //   });

  // var query_texts = svg.append("g").selectAll("text").data(dataset).enter().append("text")
  //   .text(function (d) {
  //     var p = Math.floor(Math.log(d.q) / Math.log(10));
  //     if (p==0) return 'Q = 1';
  //     if (p==1) return 'Q = 10';
  //     return 'Q = 10^' + p;
  //   })
  //   .attr({
  //     x: function(d, i) { return i * barWidth + (barWidth - barPadding) / 2; },
  //     y: function(d) { return h - hgap + 15; },
  //     'text-anchor': 'middle',
  //     'font-family': 'Arial',
  //     'font-size': '14px',
  //     // 'font-weight': 'bold',
  //   });

  // svg.append("g").append("text")
  //   .text("N = 10^8 (initially unordered), Q = Number of random point queries")
  //   .attr({
  //     x: function(d) { return w / 2; },
  //     y: function(d) { return h - 10; },
  //     'text-anchor': 'middle',
  //     'font-family': 'Arial',
  //     'font-size': '14px',
  //     'font-weight': 'bold',
  //   });
}


drawBarChartH(dataset, 1000, 400, 50, 15);
// drawBarChartV(dataset, 500, 500, 50, 15);
// var padding = 50;

// var circles = svg.selectAll("circle")
//   .data(dataset)
//   .enter()
//   .append("circle")
//   ;

// circles
//   .attr("cx", function(d, i) {
//     return (i * 50) + 25;
//   })
//   .attr("cy", h / 2)
//   .attr("r", function(d) {
//     return d.sort_time;
//   })
//   ;



// var xScale = d3.scale.linear()
//   .domain([100, 500])
//   .range([padding, w - padding])
//   ;

// var xAxis = d3.svg.axis()
//   .scale(xScale)
//   .orient("bottom")
//   .ticks(4)
//   ;

// svg.append("g")
//   .attr("class", "axis")
//   .attr("transform", "translate(0," + (h - padding) + ")")
//   .call(xAxis);


// var yScale = d3.scale.linear()
//   .domain([0, d3.max(dataset, function(d) { return d.sort_query_time; })])
//   .range([h - padding, padding])
//   ;

// var yAxis = d3.svg.axis()
//   .scale(yScale)
//   .orient("left")
//   .ticks(5);

// svg.append("g")
//   .attr("class", "axis")
//   .attr("transform", "translate(" + padding + ",0)")
//   .call(yAxis);




// d3.select("body")
//   .selectAll("div")
//   .data(dataset)
//   .enter()
//   .append("div")
//   .attr('class', 'bar_sort')
//   .style('height', function (d) {
//     return (d.sort_time + d.sort_query_time) * 2 + 'px';
//   })
//   .append("div")
//   .attr('class', 'bar_comb')
//   .style('height', function (d) {
//     return (d.comb_insert_time + d.comb_query_time) * 2 + 'px';
//   })
  // .text(function (d) {
  //   return d.q;
  // })
  // .style("color", function(d) {
  //   if (d.q > 1e4) return "red";
  //   return "black";
  // })
  ;
