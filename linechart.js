function renderLinecharts() {
  var charts = document.getElementsByClassName('chart');
  Array.prototype.forEach.call(charts, parse_div);
}

function parse_div(div, ith_div) {
  eval('var opts = ' + div.innerHTML);
  div.innerHTML = '<center>Figure ' + (ith_div + 1) + '</center>';
  linechart(div, opts);
}

function linechart(div, opts) {
  opts.width = opts.width || 400;
  opts.height = opts.height || 210;
  opts.fontSize = opts.fontSize || 14;
  opts.margin = opts.margin || { top: 20, right: 40, bottom: 45, left: 60 };

  function noFormat() { return ""; };
  noFormat.replace = noFormat;

  var superscript = "⁰¹²³⁴⁵⁶⁷⁸⁹";
  function formatPower(d) { return (d + "").split("").map(function(c) { return superscript[c]; }).join(""); };
  function ticksPow10(d) { return 10 + formatPower(Math.round(Math.log(d) / Math.LN10)); }

  var algos = opts.data = group(filter(opts.filters), opts.group_by);
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

  if (opts.xAxis.format == 'pow10') opts.xAxis.format = ticksPow10;
  if (opts.yAxis.format == 'pow10') opts.yAxis.format = ticksPow10;
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
        "d": "M" + (attr.x - 20) + " " + attr.y + " L" + (attr.x - 1) + " " + attr.y,
        "stroke": color(key),
      });
      legendG.append("path").attr({
        transform: "translate(" + (attr.x - 10) + ", " + attr.y + ")",
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
  xAxisG.selectAll("path").attr({"fill": "none", "stroke":"black"});
  xAxisG.selectAll("line").attr({"fill": "none", "stroke":"black"});

  if (opts.xAxis.scale == 'log')
    toPower10(xAxisG, opts.fontSize);

  xAxisG.append("text")
      .attr("transform", "translate(" + (width / 2) + ", 0)")
      .attr("y", opts.xAxis.margin)
      .style({"text-anchor": "middle", "font-size": opts.fontSize })
      .text(opts.xAxis.label);
  var xTopAxisG = svg.append("g")
      .attr("class", "x axis top")
      .attr("transform", "translate(0,0)")
      .style({"shape-rendering": "crispEdges", "font-size": opts.fontSize})
      .call(xAxisTop);
  xTopAxisG.selectAll("text").text(null);
  xTopAxisG.selectAll("path").attr({"fill": "none", "stroke":"black"});
  xTopAxisG.selectAll("line").attr({"fill": "none", "stroke":"black"});



  var yAxisG = svg.append("g")
      .attr("class", "y axis")
      .style({"shape-rendering": "crispEdges", "font-size": opts.fontSize})
      .call(yAxis);
  yAxisG.selectAll("text").attr("x", -5);
  yAxisG.selectAll("path").attr({"fill": "none", "stroke":"black"});
  yAxisG.selectAll("line").attr({"fill": "none", "stroke":"black"});

  if (opts.yAxis.scale == 'log' && opts.yAxis.format == 'pow10')
    toPower10(yAxisG, opts.fontSize);

  yAxisG.append("text")
      .attr("transform", "rotate(-90) translate(-" + 1*(height/2) + ", 0)")
      .attr("y", -opts.yAxis.margin)
      .style({"text-anchor": "middle", "font-size": opts.fontSize })
      .text(opts.yAxis.label);
  var yRightAxisG = svg.append("g")
      .attr("class", "y axis right")
      .attr("transform", "translate(" + width + ",0)")
      .style({"shape-rendering": "crispEdges", "font-size": opts.fontSize})
      .call(yAxisRight);
  yRightAxisG.selectAll("text").text(null);
  yRightAxisG.selectAll("path").attr({"fill": "none", "stroke":"black"});
  yRightAxisG.selectAll("line").attr({"fill": "none", "stroke":"black"});

  var lines = svg.selectAll(".lines").data(algo_arr).enter().append("g").attr("class", "lines");
  lines.append("path")
      .attr("class", "line")
      .attr("d", function(d) { return line(d.values); })
      .style("stroke", function(d) { return color(d.key); });

  svg.selectAll(".line").attr({
    "fill": "none",
    "stroke": "steelblue",
    "stroke-width": "1.5px",
  });

  // lines.append("text")
  //     .datum(function(d) { return { key: d.key, value: d.values[d.values.length - 1] }; })
  //     .attr("transform", function(d) { return "translate(" + x(d.value[opts.xAxis.attr]) + "," + y(d.value[opts.yAxis.attr]) + ")"; })
  //     .attr("x", 3)
  //     .attr("y", 5)
  //     .style({ "font-size": opts.fontSize })
  //     .text(function(d) { return algo_name[d.key].name; });
}

function toPower10(g, fontSize) {
  var texts = g.selectAll(".tick text").text(null);
  texts.filter(isOneOrTen).text(function (d) { return d; })

  texts.filter(isPowerOfTen)
      .text(10)
    .append("tspan")
      .attr("dy", "-.5em")
      .style({"shape-rendering": "crispEdges", "font-size": fontSize - 3})
      .text(function(d) { return Math.round(Math.log(d) / Math.LN10); });
}

function isPowerOfTen(d) {
  if (isOneOrTen(d)) return false;
  return d / Math.pow(10, Math.ceil(Math.log(d) / Math.LN10 - 1e-12)) === 1;
}

function isOneOrTen(d) {
  return (d == 1 || d == 10);
}




data.forEach(function (d) {
  d.Q = parseInt(d.Q);
  d.insert_time = parseFloat(d.insert_time);
  d.query_time = parseFloat(d.query_time) + 1e-9;
  d.update_time = parseFloat(d.update_time) + 1e-9;
  d.total_time = d.insert_time + d.query_time + d.update_time;
  d.avg_time = d.total_time / d.Q;
  d.qps = d.Q / d.total_time;
});

var algo_name = {
  comb:            { name: "COMB", symbol: "diamond", color: "orange" },
  combtr:          { name: "COMB-TR", symbol: "diamond", color: "green" },
  combtr2:         { name: "COMB-TR2", symbol: "diamond", color: "black" },
  comb_art:        { name: "COMB-ART", symbol: "circle", color: "magenta" },
  comb_art_ns:     { name: "COMB-ART-NS", symbol: "circle", color: "lime" },
  comb_art_1:     { name: "COMB-ART-1", symbol: "cross", color: "brown" },
  combt8192:       { name: "COMB-Tree8192", symbol: "square", color: "red" },
  combt2048:       { name: "COMB-Tree2048", symbol: "square", color: "green" },
  combt512:        { name: "COMB-Tree512", symbol: "square", color: "blue" },
  comb800:         { name: "COMB800", symbol: "diamond", color: "red" },
  comb1600:        { name: "COMB1600", symbol: "diamond", color: "blue" },
  comb3200:        { name: "COMB3200", symbol: "diamond", color: "green" },
  comb_count:      { name: "COMB", symbol: "diamond", color:"orange" },
  art:             { name: "ART", symbol: "triangle-up", color: "lime", },
  art_best:        { name: "ARTC", symbol: "cross", color: "red", },
  art_best_eager:  { name: "ARTB", symbol: "cross", color: "green", },
  sort:            { name: "Sort", symbol: "triangle-down", color: "blue", },
  crack:           { name: "Crack", symbol: "circle", color: "red", },
  crack_count:     { name: "Crack", symbol: "circle", color: "red", },
  mdd1r:           { name: "Scrack", symbol: "triangle-down", color: "red", },
  ctree_32_64:     { name: "CT64", symbol: "square", color: "brown", },
  ctree_32_1024:   { name: "CT1024", symbol: "triangle-up", color: "magenta", },
  ctree_32_4096:   { name: "CT4096", symbol: "diamond", color: "blue", },
  ctree_eager:     { name: "BTree", symbol: "square", color: "black", },
  btree_google:    { name: "BTree", symbol: "square", color: "black", },
};

function filter(filters) {
  var ret = data;
  filters.forEach(function (filter) {
    ret = ret.filter(function (d) {
      return filter.values.indexOf(d[filter.attr]) != -1;
    });
  });
  return ret;
}

function group(arr, by) {
  var keys = {}, groups = [];
  arr.forEach(function (d) {
    if (!keys[d[by]]) keys[d[by]] = [];
    keys[d[by]].push(d);
  });
  return keys;
  // for (var i in keys) if (keys.hasOwnProperty(i)) {
  //   var a = keys[i];
  //   keys[i] = [];
  //   groups.push({ group: i, label: a.label, x: a.x, y: a.y });
  // }
}

function expTime(t) {
  function f(x) {
    var v = d3.format(".2f")(x) + "";
    while (v.length > 0 && v[v.length - 1] == '0') v = v.substring(0, v.length - 1);
    if (v.length > 0 && v[v.length - 1] == '.') v = v.substring(0, v.length - 1);
    console.log(v);
    return v;
  }
  if (t >= 3600) return f(t / 3600) + 'h';
  if (t >= 60) return f(t / 60) + 'm';
  if (t >= 1) return f(t) + 's';
  if (t >= 1e-3) return f(t*1e3) + 'ms';
  if (t >= 1e-6) return f(t*1e6) + 'µs';
  return "?";
};
















function barchart(id, xcap, ylabel, data, update_w, Q, algo_name) {
  var algos = {};
  for (var i in algo_name) if (algo_name.hasOwnProperty(i))
    algos[i] = [];

  data.forEach(function (d) {
    if (algos[d.algorithm] && d.Q == Q && d.update_workload == update_w && d.query_workload == 'Random') {
      d.Q = parseInt(d.Q);
      d.insert_time = parseFloat(d.insert_time);
      d.query_time = parseFloat(d.query_time) + 1e-9;
      d.update_time = parseFloat(d.update_time) + 1e-9;
      d.total_time = d.insert_time + d.query_time;
      if (update_w == 'LFHV') d.total_time += d.update_time;
      algos[d.algorithm].push(d);
    }
  });

  // console.log(algos);

  var svg = d3.select(id);
  if (!svg[0][0]) return;
  var margin = { top: 20, right: 0, bottom: 60, left: 60 };
  if (!ylabel) margin.left = 35;
  var width = svg.attr('width') - (margin.left + margin.right);
  var height = svg.attr('height') - (margin.top + margin.bottom);

  var x = d3.scale.ordinal().rangeRoundBands([0, width], .1);
  var y = d3.scale.linear().rangeRound([height, 0]);

  var color = d3.scale.ordinal().range(["#000", "#888", "#bbb"]);
  var xAxis = d3.svg.axis().scale(x).orient("bottom");
  var yAxis = d3.svg.axis().scale(y).orient("left").tickFormat(d3.format("1s")).ticks(5);

  svg = svg.append("g").attr("transform", "translate(" + margin.left + "," + margin.top + ")");
  color.domain(['insert', 'update', 'query']);

  var arr = d3.keys(algos).map(function(name) {
    var r = algos[name][0];
    var values = [
      { name:'insert', y0: 0, y1: r.insert_time },
      { name:'query', y0: r.insert_time, y1: r.total_time },
    ];
    if (update_w == 'LFHV') {
      values[1] = { name:'update', y0: r.insert_time, y1: r.insert_time + r.update_time };
      values[2] = { name:'query', y0: r.insert_time + r.update_time, y1: r.total_time };
    }
    // console.log(name + ' ' + JSON.stringify(values));

    return {
      name: name,
      values: values
    };
  });

  x.domain(arr.map(function(d) { return algo_name[d.name]; }));
  y.domain([0, d3.max(arr, function(d) { return d.values[d.values.length - 1].y1; })]);

  var xx = svg.append("g")
      .attr("class", "x axis")
      .attr("transform", "translate(0," + height + ")")
      .call(xAxis);
  xx.selectAll("path").attr({"fill": "none", "stroke":"black"});
  xx.selectAll("line").attr({"fill": "none", "stroke":"black"});

  xx.selectAll("text").attr({
    "font-size": 12,
    "transform": "rotate(-25) translate(-7,0)",
    "text-anchor": "end"
  });
  xx.append("text")
    .attr("transform", "translate(" + width/2 + ",0)")
    .attr("y", 45)
    // .attr("dy", ".71em")
    .style("text-anchor", "middle")
    .text(xcap);

  var yy = svg.append("g")
      .attr("class", "y axis")
      .call(yAxis);
  yy.selectAll("path").attr({"fill": "none", "stroke":"black"});
  yy.selectAll("line").attr({"fill": "none", "stroke":"black"});
  if (ylabel) yy.append("text")
      .attr("transform", "rotate(-90) translate(-" + height/2 + ",0)")
      .attr("y", -45)
      // .attr("dy", ".71em")
      .style("text-anchor", "middle")
      .text(ylabel);

  var state = svg.selectAll(".state")
      .data(arr)
    .enter().append("g")
      .attr("class", "g")
      .attr("transform", function(d) { return "translate(" + x(d.name) + ",0)"; });

  state.selectAll("rect")
      .data(function(d) { return d.values; })
    .enter().append("rect")
      .attr("width", x.rangeBand())
      .attr("y", function(d) { return y(d.y1); })
      .attr("height", function(d) { return y(d.y0) - y(d.y1); })
      .style("fill", function(d) { return color(d.name); });

  // var legend = svg.selectAll(".legend")
  //     .data(['insert', 'query'])
  //   .enter().append("g")
  //     .attr("class", "legend")
  //     .attr("transform", function(d, i) { return "translate(0," + i * 20 + ")"; });

  // legend.append("rect")
  //     .attr("x", width - 18)
  //     .attr("width", 18)
  //     .attr("height", 18)
  //     .style("fill", color);

  // legend.append("text")
  //     .attr("x", width - 24)
  //     .attr("y", 9)
  //     .attr("dy", ".35em")
  //     .style("text-anchor", "end")
  //     .text(function(d) { return d; });
}
