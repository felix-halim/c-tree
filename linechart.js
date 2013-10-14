function linechart(div) {
  eval('var opts = ' + div.innerText);
  div.innerHTML = '';

  opts.width = opts.width || 450;
  opts.height = opts.height || 200;
  opts.fontSize = opts.fontSize || 14;
  opts.margin = opts.margin || { top: 20, right: 80, bottom: 30, left: 80 };

  function noFormat() { return ""; };
  noFormat.replace = noFormat;

  var superscript = "⁰¹²³⁴⁵⁶⁷⁸⁹";
  function formatPower(d) { return (d + "").split("").map(function(c) { return superscript[c]; }).join(""); };
  function ticksPow10(d) { return 10 + formatPower(Math.round(Math.log(d) / Math.LN10)); }

  var algos = opts.data;
  if (opts.base) {
    var base = algos[opts.base];
    if (!base) alert('base not found: ' + opts.base);
    for (var algo in algos) if (algos.hasOwnProperty(algo)) {
      var arr = algos[algo];
      for (var i = 0; i < base.length; i++) {
        arr[i].ratio = arr[i][opts.yAxis.attr] / base[i][opts.yAxis.attr];
      }
    }
    // delete algos[opts.base];
    for (var i = 0; i < base.length; i++) {
      base[i].ratio = 1.0;
    }
    opts.yAxis.attr = 'ratio';
  }

  if (opts.xAxis.domain) {
    for (var algo in algos) if (algos.hasOwnProperty(algo)) {
      var arr = algos[algo].filter(function (a) {
        return opts.xAxis.domain[0] <= a[opts.xAxis.attr] && a[opts.xAxis.attr] <= opts.xAxis.domain[1];
      });
      algos[algo] = arr;
    }
  }
  var width = opts.width - (opts.margin.left + opts.margin.right);
  var height = opts.height - (opts.margin.top + opts.margin.bottom);

  var x = ((opts.xAxis.scale == 'log') ? d3.scale.log() : d3.scale.linear()).range([0, width]);
  var y = ((opts.yAxis.scale == 'log') ? d3.scale.log() : d3.scale.linear()).range([height, 0]);

  var xAxis = d3.svg.axis().scale(x).orient("bottom").ticks(20, opts.xAxis.format).tickValues(opts.xAxis.ticks).tickSize(-5);
  var xAxisTop = d3.svg.axis().scale(x).orient("top").ticks(20, noFormat).tickValues(opts.xAxis.ticks).tickSize(-5);
  var yFormat = (opts.yAxis.scale == 'log') ? opts.yAxis.format : noFormat;
  var yAxis = d3.svg.axis().scale(y).orient("left").ticks(10, yFormat).tickValues(opts.yAxis.ticks).tickSize(-5);
  var yAxisRight = d3.svg.axis().scale(y).orient("right").ticks(10, noFormat).tickValues(opts.yAxis.ticks).tickSize(-5);
  
  function fx(d) { return x(d[opts.xAxis.attr]); }
  function fy(d) { return y(d[opts.yAxis.attr]); }
  var line = d3.svg.line().interpolate("linear").x(fx).y(fy);

  // var color = d3.scale.category10();
  var color = d3.scale.ordinal();
  color.domain(d3.keys(algos));
  color.range(color.domain().map(function (key) { return algo_name[key].color; }));

  var algo_arr = color.domain().map(function(key) { return { key: key, values: algos[key] }; });
  function extreme(f, arr, attr) { return f(arr, function (a) { return f(a.values, function (d) { return d[attr]; }); }); }
  // x.domain(d3.extent(data, function(d) { return d.Q; }));
  x.domain(opts.xAxis.domain || [ extreme(d3.min, algo_arr, opts.xAxis.attr), extreme(d3.max, algo_arr, opts.xAxis.attr) ]);
  y.domain(opts.yAxis.domain || [ extreme(d3.min, algo_arr, opts.yAxis.attr), extreme(d3.max, algo_arr, opts.yAxis.attr) ]);

  var svg = d3.select(div).append('svg');
  svg.attr('width', opts.width);
  svg.attr('height', opts.height);
  div.style.width = opts.width+'px';
  div.style.height = opts.height+'px';
  svg = svg.append("g").attr("transform", "translate(" + opts.margin.left + "," + opts.margin.top + ")");

  function symbol(algo) { return d3.svg.symbol().size(30).type(algo_name[algo].symbol); }

  var arr = [];
  color.domain().map(function(key) { arr = arr.concat(algos[key]); });
  svg.selectAll("path.dot")
     .data(arr)
     .enter().append("path")
     .attr("class", "dot")
     .attr("transform", function (d) { return "translate(" + fx(d) + ", " + fy(d) + ")"; })
     .attr("d", function (d) { return symbol(d.algorithm)(); })
     .attr("fill", function (d) { return color(d.algorithm); })
     .attr("stroke", function (d) { return color(d.algorithm); });

  if (opts.legend) {
    var legendG = svg.append("g");
    d3.keys(opts.legend).forEach(function (key) {
      var attr = opts.legend[key];
      attr.fill = color(key);
      attr.style = "font-size:14px";
      attr["alignment-baseline"] = 'middle';
      legendG.append("text").attr(attr).text(algo_name[key].name);
      legendG.append("path").attr({
        "d": "M" + (attr.x - 25) + " " + attr.y + " L" + (attr.x - 5) + " " + attr.y,
        "stroke": color(key),
      });
      legendG.append("path").attr({
        transform: "translate(" + (attr.x - 15) + ", " + attr.y + ")",
        fill: color(key),
        d: symbol(key),
      });
    });
  }

  opts.fontSize = opts.fontSize || "15px";

  var xAxisG = svg.append("g")
      .attr("class", "x axis")
      .attr("transform", "translate(0," + height + ")")
      .style({"shape-rendering": "crispEdges", "font-size": opts.fontSize})
      .call(xAxis);
  xAxisG.selectAll("text").attr("y", 7);
  if (opts.xAxis.scale == 'log')
    toPower10(xAxisG, opts.fontSize);

  xAxisG.append("text")
      .attr("transform", "translate(" + (width / 2) + ", 0)")
      .attr("y", opts.xAxis.margin)
      .style({"text-anchor": "middle", "font-size": opts.fontSize })
      .text(opts.xAxis.label);
  svg.append("g")
      .attr("class", "x axis top")
      .attr("transform", "translate(0,0)")
      .style({"shape-rendering": "crispEdges", "font-size": opts.fontSize})
      .call(xAxisTop)
      .selectAll("text")
      .text(null);



  var yAxisG = svg.append("g")
      .attr("class", "y axis")
      .style({"shape-rendering": "crispEdges", "font-size": opts.fontSize})
      .call(yAxis);
  yAxisG.selectAll("text").attr("x", -5);
  if (opts.yAxis.scale == 'log')
    toPower10(yAxisG, opts.fontSize);

  yAxisG.append("text")
      .attr("transform", "rotate(-90) translate(-" + 1*(height/2) + ", 0)")
      .attr("y", -opts.yAxis.margin)
      .style({"text-anchor": "middle", "font-size": opts.fontSize })
      .text(opts.yAxis.label);
  svg.append("g")
      .attr("class", "y axis right")
      .attr("transform", "translate(" + width + ",0)")
      .style({"shape-rendering": "crispEdges", "font-size": opts.fontSize})
      .call(yAxisRight)
      .selectAll("text")
      .text(null);

  var lines = svg.selectAll(".lines").data(algo_arr).enter().append("g").attr("class", "lines");
  lines.append("path")
      .attr("class", "line")
      .attr("d", function(d) { return line(d.values); })
      .style("stroke", function(d) { return color(d.key); });

  // lines.append("text")
  //     .datum(function(d) { return { key: d.key, value: d.values[d.values.length - 1] }; })
  //     .attr("transform", function(d) { return "translate(" + x(d.value[opts.xAxis.attr]) + "," + y(d.value[opts.yAxis.attr]) + ")"; })
  //     .attr("x", 3)
  //     .attr("y", 5)
  //     .style({ "font-size": opts.fontSize })
  //     .text(function(d) { return algo_name[d.key].name; });
}

function toPower10(g, fontSize) {
  g.selectAll(".tick text")
      .text(null)
    .filter(isPowerOfTen)
      .text(10)
    .append("tspan")
      .attr("dy", "-.5em")
      .style({"shape-rendering": "crispEdges", "font-size": fontSize - 3})
      .text(function(d) { return Math.round(Math.log(d) / Math.LN10); });
}

function isPowerOfTen(d) {
  return d / Math.pow(10, Math.ceil(Math.log(d) / Math.LN10 - 1e-12)) === 1;
}
