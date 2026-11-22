// Diffusor documentation site: theme switch, heading anchors, the
// "On this page" list, scrollable tables and copy buttons on code blocks.
(function () {
  "use strict";

  var root = document.documentElement;
  var prose = document.querySelector(".prose");

  // ----- light / dark switch (remembers the choice when storage works)
  var toggle = document.querySelector(".theme-toggle");
  if (toggle) {
    toggle.addEventListener("click", function () {
      var current = root.dataset.theme ||
        (window.matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light");
      var next = current === "dark" ? "light" : "dark";
      root.dataset.theme = next;
      try { localStorage.setItem("diffusor-theme", next); } catch (e) {}
    });
  }

  if (!prose) return;

  // ----- tables scroll sideways inside a frame instead of widening the page
  prose.querySelectorAll("table").forEach(function (table) {
    var wrap = document.createElement("div");
    wrap.className = "table-wrap";
    table.parentNode.insertBefore(wrap, table);
    wrap.appendChild(table);
  });

  // ----- copy buttons on code blocks
  if (navigator.clipboard) {
    prose.querySelectorAll("div.highlighter-rouge").forEach(function (block) {
      var code = block.querySelector("pre");
      if (!code) return;
      var button = document.createElement("button");
      button.type = "button";
      button.className = "copy-btn";
      button.textContent = "Copy";
      button.addEventListener("click", function () {
        navigator.clipboard.writeText(code.innerText.replace(/\n$/, "")).then(function () {
          button.textContent = "Copied";
          setTimeout(function () { button.textContent = "Copy"; }, 1500);
        });
      });
      block.appendChild(button);
    });
  }

  // ----- heading anchors and the "On this page" list
  var headings = Array.prototype.slice.call(prose.querySelectorAll("h2[id], h3[id]"));
  headings.forEach(function (h) {
    var a = document.createElement("a");
    a.className = "anchor";
    a.href = "#" + h.id;
    a.setAttribute("aria-label", "Link to this section");
    a.textContent = "#";
    h.insertBefore(a, h.firstChild);
  });

  var toc = document.querySelector(".toc");
  var list = toc && toc.querySelector(".toc-list");
  if (!list || headings.length < 2) return;

  var links = {};
  headings.forEach(function (h) {
    var li = document.createElement("li");
    li.className = "lvl-" + h.tagName.charAt(1);
    var a = document.createElement("a");
    a.href = "#" + h.id;
    a.textContent = h.textContent.replace(/^#/, "").trim();
    li.appendChild(a);
    list.appendChild(li);
    links[h.id] = a;
  });
  toc.hidden = false;

  // highlight the last heading that has scrolled above the top third
  var active = null;
  function update() {
    var line = window.innerHeight * 0.3;
    var current = headings[0];
    for (var i = 0; i < headings.length; i++) {
      if (headings[i].getBoundingClientRect().top <= line) current = headings[i];
      else break;
    }
    var link = links[current.id];
    if (link === active) return;
    if (active) active.classList.remove("active");
    link.classList.add("active");
    active = link;
    var box = toc.getBoundingClientRect();
    var r = link.getBoundingClientRect();
    if (r.top < box.top || r.bottom > box.bottom) {
      toc.scrollTop += r.top - box.top - box.height / 3;
    }
  }
  var ticking = false;
  window.addEventListener("scroll", function () {
    if (ticking) return;
    ticking = true;
    requestAnimationFrame(function () { ticking = false; update(); });
  }, { passive: true });
  update();
})();
