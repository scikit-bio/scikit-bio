// originally taken from scikit-learn's Sphinx theme
$(document).ready(function() {
    /* Add a [>>>] button on the top-right corner of code samples to hide
     * the >>> and ... prompts and the output and thus make the code
     * copyable.
     * Note: This JS snippet was taken from the official python.org
     * documentation site.*/
    var div = $('.highlight-python .highlight,' +
                '.highlight-python3 .highlight,' +
                '.highlight-pycon .highlight')
    var pre = div.find('pre');

    // get the styles from the current theme
    pre.parent().parent().css('position', 'relative');
    var hide_text = 'Hide the prompts and output';
    var show_text = 'Show the prompts and output';
    var border_width = pre.css('border-top-width');
    var border_style = pre.css('border-top-style');
    var border_color = pre.css('border-top-color');
    var button_styles = {
        'cursor':'pointer', 'position': 'absolute', 'top': '0', 'right': '0',
        'border-color': border_color, 'border-style': border_style,
        'border-width': border_width, 'color': border_color, 'text-size': '75%',
        'font-family': 'monospace', 'padding-left': '0.2em', 'padding-right': '0.2em',
        'display': 'inline'
    }

    // create and add the button to all the code blocks that contain >>>
    div.each(function(index) {
        var jthis = $(this);
        if (jthis.find('.gp').length > 0) {
            var button = $('<span class="copybutton">&gt;&gt;&gt;</span>');
            button.css(button_styles)
            button.attr('title', hide_text);
            jthis.prepend(button);

            var show_output = false;
            button.bind('click', function() {
                if (show_output) {
                    var button = $(this);
                    button.parent().find('.go, .gp, .gt').show();
                    button.next('pre').find('.gt').nextUntil('.gp, .go').css('visibility', 'visible');
                    button.css('text-decoration', 'none');
                    button.attr('title', hide_text);
                    show_output = false;
                } else {
                    var button = $(this);
                    button.parent().find('.go, .gp, .gt').hide();
                    button.next('pre').find('.gt').nextUntil('.gp, .go').css('visibility', 'hidden');
                    button.css('text-decoration', 'line-through');
                    button.attr('title', show_text);
                    show_output = true;
                }
            });
        }
        // tracebacks (.gt) contain bare text elements that need to be
        // wrapped in a span to work with .nextUntil() (see later)
        jthis.find('pre:has(.gt)').contents().filter(function() {
            return ((this.nodeType == 3) && (this.data.trim().length > 0));
        }).wrap('<span>');
    });
});
