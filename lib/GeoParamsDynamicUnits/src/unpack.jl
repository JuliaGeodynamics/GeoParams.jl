using UnPack

export @unpack_val, @unpack_units

"""
    @unpack_val a, b = object

Unpack the numeric values stored by `GeoUnit` fields.
"""
macro unpack_val(args)
    args.head == :(=) || error("Expression needs to be of form `a, b = c`")
    items, source = args.args
    items = items isa Symbol ? [items] : items.args
    source_instance = gensym()

    assignments = [
        :($key = $UnPack.unpack($source_instance, Val{$(QuoteNode(key))}()).val) for
            key in items
    ]

    return esc(
        quote
            local $source_instance = $source
            $(Expr(:block, assignments...))
            $source_instance
        end
    )
end

"""
    @unpack_units a, b = object

Unpack the dimensional values stored by `GeoUnit` fields.
"""
macro unpack_units(args)
    args.head == :(=) || error("Expression needs to be of form `a, b = c`")
    items, source = args.args
    items = items isa Symbol ? [items] : items.args
    source_instance = gensym()

    assignments = [
        quote
                $key =
                $UnPack.unpack($source_instance, Val{$(QuoteNode(key))}()).val .*
                $UnPack.unpack($source_instance, Val{$(QuoteNode(key))}()).unit
            end for key in items
    ]

    return esc(
        quote
            local $source_instance = $source
            $(Expr(:block, assignments...))
            $source_instance
        end
    )
end
