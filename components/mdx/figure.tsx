import Image from "next/image"

interface FigureProps {
  src: string
  caption?: string
  alt?: string
  width?: number
  height?: number
}

export function Figure({ src, caption, alt, width = 800, height = 400 }: FigureProps) {
  return (
    <figure className="my-8">
      <div className="relative overflow-hidden rounded-lg border bg-muted">
        <Image
          src={src || `/placeholder.svg?height=${height}&width=${width}&query=academic figure`}
          alt={alt || caption || ""}
          width={width}
          height={height}
          className="w-full object-cover"
          unoptimized
        />
      </div>
      {caption && (
        <figcaption className="mt-3 text-center text-sm text-muted-foreground leading-relaxed">{caption}</figcaption>
      )}
    </figure>
  )
}
