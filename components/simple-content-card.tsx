import { Card, CardContent, CardHeader } from "@/components/ui/card"
import { Badge } from "@/components/ui/badge"
import { Button } from "@/components/ui/button"
import { CalendarDays, Clock, User } from "lucide-react"
import { format } from "date-fns"
import { zhCN } from "date-fns/locale"
import Link from "next/link"
import Image from "next/image"
import type { PostMetadata } from "@/lib/simple-content-loader"

interface SimplePostCardProps {
  post: PostMetadata
}

export function SimplePostCard({ post }: SimplePostCardProps) {
  return (
    <Card className="transition-all hover:shadow-md">
      {post.coverImage && (
        <div className="aspect-video overflow-hidden rounded-t-lg">
          <Image
            src={post.coverImage || "/placeholder.svg"}
            alt={post.title}
            width={400}
            height={225}
            className="w-full h-full object-cover transition-transform hover:scale-105"
            unoptimized
          />
        </div>
      )}
      <CardHeader className="pb-3">
        <div className="flex flex-wrap gap-2 mb-3">
          {post.tags.slice(0, 3).map((tag) => (
            <Badge key={tag} variant="secondary" className="text-xs">
              {tag}
            </Badge>
          ))}
          {post.tags.length > 3 && (
            <Badge variant="outline" className="text-xs">
              +{post.tags.length - 3}
            </Badge>
          )}
        </div>
        <h3 className="text-lg font-semibold text-balance leading-tight mb-2">
          <Link href={`/${post.slug}`} className="hover:text-primary transition-colors">
            {post.title}
          </Link>
        </h3>
      </CardHeader>
      <CardContent>
        <p className="text-muted-foreground text-sm text-pretty mb-4 line-clamp-3">{post.excerpt}</p>
        <div className="flex items-center justify-between">
          <div className="flex items-center gap-4 text-xs text-muted-foreground">
            <div className="flex items-center gap-1">
              <User className="h-3 w-3" />
              <span>{post.authors[0]}</span>
              {post.authors.length > 1 && <span>等</span>}
            </div>
            <div className="flex items-center gap-1">
              <CalendarDays className="h-3 w-3" />
              <span>{format(new Date(post.publishedAt), "MM/dd", { locale: zhCN })}</span>
            </div>
          </div>
          <Button asChild variant="ghost" size="sm">
            <Link href={`/${post.slug}`}>阅读全文</Link>
          </Button>
        </div>
      </CardContent>
    </Card>
  )
}