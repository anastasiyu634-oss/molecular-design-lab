"""Best-effort Chrome/Edge page-translation resilience for Streamlit.

This module does not translate the UI. It injects an experimental client-side
DOM resilience shim into Streamlit's document, so browser translation can
continue rather than corrupting React's references to text nodes.

The JavaScript is loaded only from a pinned public third-party version.
Browser must be online and allow the CDN. See README_TRANSLATION.txt.
"""
from __future__ import annotations

import inspect
import streamlit as st

# MIT open-source dependency: https://github.com/alexspeller/translation-resilience
# Version pinned to avoid silent changes. JS is external, not shipped or edited.
_TRANSLATION_GUARD = r"""
<script>
(function () {
  // The original app renders Russian text. Streamlit's index.html declares
  // English by default, which can trigger unnecessary translation prompts.
  var root = document.documentElement;
  if (!root.classList.contains('translated-ltr') &&
      !root.classList.contains('translated-rtl') &&
      !window.__mdlabLocaleWasSet) {
    root.lang = 'ru';
    window.__mdlabLocaleWasSet = true;
  }

  // Streamlit reruns the Python app on button presses, but this DOM guard
  // should be installed only once in the current browser document.
  if (window.__mdlabTranslationGuardState === 'active' ||
      window.__mdlabTranslationGuardState === 'loading') {
    return;
  }
  window.__mdlabTranslationGuardState = 'loading';
  import('https://cdn.jsdelivr.net/npm/translation-resilience@0.3.0/+esm')
    .then(function (module) {
      if (typeof module.installTranslationResilience !== 'function') {
        throw new Error('Missing installTranslationResilience export');
      }
      // eager supports Chrome Translate already active when Streamlit starts.
      window.__mdlabTranslationGuardUninstall =
        module.installTranslationResilience({eager: true});
      window.__mdlabTranslationGuardState = 'active';
    })
    .catch(function (error) {
      window.__mdlabTranslationGuardState = 'error';
      console.warn('[Molecular Design Lab] Translation guard unavailable:', error);
    });
})();
</script>
"""


def install_browser_translation_guard() -> None:
    """Install before UI widgets; do nothing on older Streamlit versions."""
    try:
        accepts_js = 'unsafe_allow_javascript' in inspect.signature(st.html).parameters
    except (AttributeError, TypeError, ValueError):
        accepts_js = False
    if accepts_js:
        st.html(_TRANSLATION_GUARD, unsafe_allow_javascript=True)
    else:
        st.warning('Для устойчивости к переводу Chrome требуется Streamlit 1.52 или новее. '
                   'Обновите библиотеку в окружении сайта.')
