/* Interactive charts, kept light (rules in the project notes):
   - a <div class="sg-plot" data-src="....json"> holds a static picture (the fallback);
   - when it scrolls into view, the basic Plotly build is loaded once and the picture is
     replaced by the interactive chart drawn from the JSON;
   - charts follow the site's light/dark switch. */
(function () {
  "use strict";
  var LIB = "https://cdn.jsdelivr.net/npm/plotly.js-basic-dist-min@3.0.1/plotly-basic.min.js";
  var libPromise = null;

  function loadLib() {
    if (window.Plotly) return Promise.resolve(window.Plotly);
    if (!libPromise) {
      libPromise = new Promise(function (resolve, reject) {
        var s = document.createElement("script");
        s.src = LIB; s.async = true;
        s.onload = function () { resolve(window.Plotly); };
        s.onerror = reject;
        document.head.appendChild(s);
      });
    }
    return libPromise;
  }

  function dark() {
    return document.body.getAttribute("data-md-color-scheme") === "slate";
  }

  function colours(layout) {
    var d = dark();
    var ink = d ? "#e8e8ec" : "#0b0b0b", ink2 = d ? "#c4c4cc" : "#52514e";
    var grid = d ? "#2c2f34" : "#e1e0d9", axis = d ? "#4a4d55" : "#c3c2b7";
    var surface = d ? "#16181b" : "#fcfcfb";
    var u = { paper_bgcolor: surface, plot_bgcolor: surface, "font.color": ink2,
              "title.font.color": ink, "legend.font.color": ink2 };
    Object.keys(layout).forEach(function (k) {
      if (/^[xy]axis\d*$/.test(k)) {
        u[k + ".gridcolor"] = grid; u[k + ".linecolor"] = axis;
        u[k + ".tickfont.color"] = ink2; u[k + ".title.font.color"] = ink2;
      }
    });
    (layout.annotations || []).forEach(function (a, i) { u["annotations[" + i + "].font.color"] = ink2; });
    return u;
  }

  function draw(box) {
    if (box.dataset.state) return;
    box.dataset.state = "loading";
    Promise.all([loadLib(), fetch(box.dataset.src).then(function (r) {
      if (!r.ok) throw new Error(r.status); return r.json();
    })]).then(function (res) {
      var Plotly = res[0], spec = res[1], layout = spec.layout || {};
      delete layout.width; layout.autosize = true;
      var target = document.createElement("div");
      target.className = "sg-plot-live";
      box.replaceChildren(target);
      return Plotly.newPlot(target, spec.data, layout,
        { responsive: true, displaylogo: false, modeBarButtonsToRemove: ["lasso2d", "select2d"] })
        .then(function () { Plotly.relayout(target, colours(target.layout)); box.dataset.state = "live"; });
    }).catch(function () { box.dataset.state = "picture"; });   /* the picture stays */
  }

  function recolour() {
    if (!window.Plotly) return;
    document.querySelectorAll(".sg-plot-live").forEach(function (el) {
      window.Plotly.relayout(el, colours(el.layout));
    });
  }

  function init() {
    var boxes = document.querySelectorAll(".sg-plot[data-src]");
    if (!boxes.length) return;
    if (!("IntersectionObserver" in window)) { boxes.forEach(draw); return; }
    var io = new IntersectionObserver(function (entries) {
      entries.forEach(function (e) { if (e.isIntersecting) { io.unobserve(e.target); draw(e.target); } });
    }, { rootMargin: "300px" });
    boxes.forEach(function (b) { io.observe(b); });
  }

  new MutationObserver(recolour).observe(document.body, { attributes: true, attributeFilter: ["data-md-color-scheme"] });
  if (window.document$ && window.document$.subscribe) { window.document$.subscribe(init); }   /* Material page loads */
  else { document.addEventListener("DOMContentLoaded", init); }
})();
