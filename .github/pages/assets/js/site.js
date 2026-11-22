// Diffusor documentation site: theme switch, the download list, heading
// anchors, the "On this page" list, scrollable tables and copy buttons on code
// blocks.
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

  // ----- the Download button opens the list of installers instead of GitHub
  var dlToggle = document.querySelector(".dl-toggle");
  var dlPanel = document.getElementById("downloads");
  if (dlToggle && dlPanel) setUpDownloads(dlToggle, dlPanel);

  function setUpDownloads(button, panel) {
    // the visitor's system, to tag its row; two Mac or two Linux files are both tagged
    var ua = navigator.userAgent || "";
    var platform = (navigator.userAgentData && navigator.userAgentData.platform) ||
      navigator.platform || "";
    var os = /Win/i.test(platform) || /Windows/i.test(ua) ? "windows"
      : /Mac/i.test(platform) && !/iPhone|iPad/i.test(ua) ? "mac"
      : /Linux|X11/i.test(platform) && !/Android/i.test(ua) ? "linux" : "";
    panel.querySelectorAll("li[data-os]").forEach(function (li) {
      if (li.dataset.os === os) li.classList.add("is-yours");
    });

    var loaded = false;
    function open(show) {
      panel.hidden = !show;
      button.setAttribute("aria-expanded", String(show));
      if (show && !loaded) { loaded = true; latestRelease(); }
      if (show) {
        // on a short screen the list runs below the fold; scroll just enough to show it
        var calm = window.matchMedia("(prefers-reduced-motion: reduce)").matches;
        panel.scrollIntoView({ block: "nearest", behavior: calm ? "auto" : "smooth" });
      }
    }
    button.addEventListener("click", function (ev) {
      ev.preventDefault();
      open(panel.hidden);
    });
    document.addEventListener("keydown", function (ev) {
      if (ev.key === "Escape" && !panel.hidden) { open(false); button.focus(); }
    });
    document.addEventListener("click", function (ev) {
      if (!panel.hidden && !panel.contains(ev.target) && !button.contains(ev.target)) open(false);
    });
    // the "How to open it" link jumps to the install notes further down the page
    panel.querySelectorAll('a[href^="#"]').forEach(function (a) {
      a.addEventListener("click", function () { open(false); });
    });

    // The links in the page point at the version the site was built with. The
    // latest release is asked for when the list first opens, so the links and
    // sizes follow a release made after the site was built. If GitHub does not
    // answer, the built links stay.
    function latestRelease() {
      var repo = panel.dataset.repo;
      if (!repo || !window.fetch) return;
      fetch("https://api.github.com/repos/" + repo + "/releases/latest",
            { headers: { Accept: "application/vnd.github+json" } })
        .then(function (r) { return r.ok ? r.json() : null; })
        .then(function (release) {
          if (!release || !release.assets) return;
          panel.querySelectorAll("a[data-suffix]").forEach(function (a) {
            var asset = release.assets.filter(function (x) {
              return x.name.slice(-a.dataset.suffix.length) === a.dataset.suffix;
            })[0];
            if (!asset) return;
            a.href = asset.browser_download_url;
            var size = a.querySelector(".dl-size");
            if (size) size.textContent = Math.round(asset.size / 1048576) + " MB";
          });
          var label = button.querySelector(".dl-label");
          var tag = String(release.tag_name || "").replace(/^v/, "");
          if (label && tag) label.textContent = "Download v" + tag;
        })
        .catch(function () {});
    }
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
