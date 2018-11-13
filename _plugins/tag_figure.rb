module Jekyll
  class FigureTag < Liquid::Block
    def initialize(tag_name, markup, tokens)
      super
      @attributes = {}

      markup.scan(Liquid::TagAttributes) do |key, value|
        @attributes[key] = value.gsub(/^'|"/, '').gsub(/'|"$/, '')
      end
    end

    def render(context)
      site = context.registers[:site]
      converter = site.find_converter_instance(::Jekyll::Converters::Markdown)

      figure_classes = @attributes['class'].to_s
      figure_main = converter.convert(super(context))
      figure_caption = converter.convert(@attributes['caption'].to_s)
      figure_topcaption = converter.convert(@attributes['topcaption'].to_s)
      # Render <figure>
      if @attributes['class']
        source = "<figure class=\"#{figure_classes}\">"
      else
        source = "<figure class=\"figure\">"
      end

      if @attributes['topcaption']
        source += "<figcaption class=\"figure__caption\">#{figure_topcaption}</figcaption>"
      end
      source += "<div class=\"figure__main\">#{figure_main}</div>"

      if @attributes['caption']
        source += "<figcaption class=\"figure__caption\" style=\"padding-top:0;\">#{figure_caption}</figcaption>"
      end
      source += "</figure>"
    end
  end
end

Liquid::Template.register_tag('figure', Jekyll::FigureTag)
