// static/js/main.js

// Import application state
import "./utils/app-state.js";

// Import form handlers
import FormResetHandler from "./handlers/form-reset-handler.js";
import FormSubmissionHandler from "./handlers/form-submission-handler.js";

// Import UI components
import DarkModeManager from "./ui/dark-mode-manager.js";
import CharacterCountManager from "./ui/character-count-manager.js";
import TabScroller from "./ui/sequence-tabs/sequence-tabs-scroller.js";
import SequenceTabsManager from "./ui/sequence-tabs/sequence-tabs-manager.js";
import TooltipManager from "./ui/tooltip-manager.js";

// Import utilities
import "./utils/utils.js";
import "./utils/validation.js";

document.addEventListener("DOMContentLoaded", () => {
    console.log("🚀 Frontend reporting for business!");
    console.log("📅 Page load time:", new Date().toISOString());
    console.log("🔧 Starting initialization process...");

    // Initialize core components
    console.log("🧩 Initializing core components");
    
    console.log("🔄 Creating SequenceTabsManager");
    const sequenceTabsManager = new SequenceTabsManager();
    window.sequenceTabsManagerInstance = sequenceTabsManager;
    console.log("✅ SequenceTabsManager created and assigned to window");

    console.log("📝 Creating FormSubmissionHandler");
    const formSubmissionHandler = new FormSubmissionHandler();
    console.log("✅ FormSubmissionHandler created");
    
    console.log("🔄 Creating FormResetHandler");
    const formResetHandler = new FormResetHandler();
    console.log("✅ FormResetHandler created");
    
    console.log("🔄 Creating TabScroller");
    new TabScroller("sequence-tabs-container");
    console.log("✅ TabScroller created");
    
    console.log("🔄 Creating CharacterCountManager");
    new CharacterCountManager();
    console.log("✅ CharacterCountManager created");
    
    console.log("🔄 Creating TooltipManager");
    new TooltipManager();
    console.log("✅ TooltipManager created");
    
    console.log("🔄 Creating DarkModeManager");
    new DarkModeManager();
    console.log("✅ DarkModeManager created");

    // Update character counters for pre-filled fields
    setTimeout(() => {
        document.querySelectorAll(".dynamic-sequence-input").forEach(input => {
            if (input.value) {
                sequenceTabsManager.updateCharCount?.(input);
            }
        });
    }, 50);

    const verboseModeCheckbox = document.getElementById("verbose_mode");
    if (verboseModeCheckbox) {
        verboseModeCheckbox.checked = true;
    }

    // Sequence count display and controls
    const numDisplay = document.getElementById("numDisplay");
    const updateNumDisplay = () => {
        if (numDisplay) {
            numDisplay.textContent = sequenceTabsManager.currentCount;
        }
    };

    document.getElementById("incrementBtn")?.addEventListener("click", () => {
        sequenceTabsManager.incrementSequenceCount();
        updateNumDisplay();
    });

    document.getElementById("decrementBtn")?.addEventListener("click", () => {
        sequenceTabsManager.decrementSequenceCount();
        updateNumDisplay();
    });

    // CENTRALIZED EVENT LISTENERS
    // ==========================
    console.log('🔍 Setting up centralized event listeners');

    // Form submission handler
    const form = document.getElementById("primer-form");
    if (form) {
        console.log('📝 Form found, attaching submit listener');
        form.addEventListener('submit', (event) => {
            console.log('🚀 FORM SUBMIT EVENT TRIGGERED');
            const submitButton = form.querySelector("[type='submit']");
            console.log('Submit button state:', {
                disabled: submitButton.disabled,
                processing: submitButton.classList.contains('processing')
            });
            
            // Prevent default form submission if already processing
            if (submitButton.disabled && submitButton.classList.contains('processing')) {
                console.log('⚠️ Form already processing, preventing submission');
                event.preventDefault();
                event.stopPropagation();
                return false;
            }
            
            // Disable button and show loading state
            console.log('⏳ Disabling submit button and showing loading state');
            formSubmissionHandler.disableSubmitButton();
            submitButton.classList.add('processing');
            submitButton.innerHTML = '<span class="spinner"></span> Processing...';
            console.log('Submit event processing complete, continuing with form submission');
        });
    } else {
        console.warn('⚠️ Form element not found in DOM');
    }

    // HTMX global event listeners
    if (typeof htmx !== 'undefined') {
        console.log('✅ HTMX is loaded, version:', htmx.version);
        
        // Central request configuration
        console.log('📡 Setting up htmx:configRequest listener');
        document.body.addEventListener('htmx:configRequest', (evt) => {
            console.log('🔄 HTMX CONFIG REQUEST EVENT FIRED', new Date().toISOString());
            
            const elt = evt.detail.elt;
            console.log('Element that triggered request:', elt);
            
            // For regular form submissions, let HTMX handle the form data naturally
            if (elt.id === 'primer-form') {
                console.log('📋 Form submission detected - using native HTMX form handling');
                // DO NOT manually set parameters - let HTMX handle the form normally
                console.log('Using HTMX default form handling');
            }
        });

        // Validation response handler
        console.log('📡 Setting up htmx:afterOnLoad listener');
        document.body.addEventListener('htmx:afterOnLoad', (evt) => {
            console.log('📥 HTMX AFTER LOAD EVENT FIRED', new Date().toISOString());
            console.log('Response element:', evt.detail.elt);
            console.log('hx-post attribute:', evt.detail.elt.getAttribute('hx-post'));
            
            if (evt.detail.elt.getAttribute('hx-post')?.includes('validate_field')) {
                console.log('✅ Field validation response received, updating submit button state');
                formSubmissionHandler.updateSubmitButtonState();
                console.log('Submit button state updated');
            } else {
                console.log('Not a validation response, no submit button update needed');
            }
            
            console.log('htmx:afterOnLoad handler complete');
        });

        // Request logging - consolidated from beforeRequest and beforeSend
        console.log('📡 Setting up htmx:beforeRequest listener');
        document.body.addEventListener('htmx:beforeRequest', (evt) => {
            console.log('📤 HTMX BEFORE REQUEST EVENT FIRED', new Date().toISOString());
            console.log('Request target URL:', evt.detail.pathInfo.requestPath);
            console.log('Request method:', evt.detail.verb);
            console.log('Request parameters:', evt.detail.parameters);
            console.log('Request headers:', evt.detail.headers);
            console.log('Triggering element:', evt.detail.elt);
            console.log('⏩ Request is about to be sent to server');
        });

        // Response logging - consolidated from afterRequest and xhr:loadend
        console.log('📡 Setting up htmx:afterRequest listener');
        document.body.addEventListener('htmx:afterRequest', (evt) => {
            console.log('📩 HTMX AFTER REQUEST EVENT FIRED', new Date().toISOString());
            console.log('Request completed to URL:', evt.detail.pathInfo.requestPath);
            console.log('Response status:', evt.detail.xhr.status);
            console.log('Response status text:', evt.detail.xhr.statusText);
            
            try {
                const responseJson = JSON.parse(evt.detail.xhr.response);
                console.log('Response JSON:', responseJson);
            } catch (e) {
                console.log('Response (not JSON):', evt.detail.xhr.response);
                console.log('Response text:', evt.detail.xhr.responseText);
            }
            
            console.log('Response headers:', evt.detail.xhr.getAllResponseHeaders());
            console.log('Request cycle complete');
        });

        // Error handling for HTMX requests
        console.log('📡 Setting up htmx:responseError listener');
        document.body.addEventListener('htmx:responseError', (evt) => {
            console.error('❌ HTMX RESPONSE ERROR EVENT FIRED', new Date().toISOString());
            console.error('Request failed to URL:', evt.detail.pathInfo?.requestPath);
            console.error('Error details:', {
                error: evt.detail.error,
                status: evt.detail.xhr.status,
                statusText: evt.detail.xhr.statusText,
                response: evt.detail.xhr.response
            });
            
            try {
                const responseJson = JSON.parse(evt.detail.xhr.response);
                console.error('Error response JSON:', responseJson);
            } catch (e) {
                console.error('Error response (not JSON):', evt.detail.xhr.response);
            }
            
            console.error('Error handling complete');
        });
        
        // Additional event listener for debugging
        console.log('📡 Setting up htmx:beforeSend listener for extra debugging');
        document.body.addEventListener('htmx:beforeSend', (evt) => {
            console.log('🔍 HTMX BEFORE SEND EVENT FIRED', new Date().toISOString());
            console.log('XHR object:', evt.detail.xhr);
            console.log('Final request parameters:', evt.detail.parameters);
            console.log('Final request headers:', {...evt.detail.headers});
        });
        
        // Additional event listener for request start
        console.log('📡 Setting up htmx:sendError listener for request failures');
        document.body.addEventListener('htmx:sendError', (evt) => {
            console.error('⛔ HTMX SEND ERROR EVENT FIRED', new Date().toISOString());
            console.error('Request failed to send:', evt.detail);
            console.error('Error:', evt.detail.error);
        });
    } else {
        console.error('⛔ HTMX is not loaded! AJAX functionality will not work');
    }
});